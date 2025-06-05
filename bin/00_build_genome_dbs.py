#!/usr/bin/env python

#$ -jc long
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID
#$ -cwd
#$ -N build_genome_dbs.py

"""
Creates a database of all 'complete' genome sequences for a TAXID, reformats
them to fasta, while capturing metadata from biosample and the genome record
"""

import sys
import gzip
import argparse
from pathlib import Path
from datetime import datetime
import json
import re

from Bio import SeqIO
from lxml import etree
from tqdm import tqdm
import pandas as pd
import requests

sys.path.append('bin')
from common import (
    get_taxa_name,
    clean_description,
    make_request,
)

NCBI_TAXID = "1423"
ENA_URI = "https://www.ebi.ac.uk/ena/"
ENVO_URI = "https://www.ebi.ac.uk/ols4/api/ontologies/envo/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F"


def search_available():
    """
    Carries out search of ENA using portal api to identify available genomes

    required parameters:
        None

    returns:
        genome_info(dict)
    """

    uri = f"{ENA_URI}portal/api/search?result=assembly&query=tax_tree({NCBI_TAXID})&fields=assembly_title,tax_id&format=json"
    result = make_request(uri)
    json_data = json.loads(result)

    # The search returns multiple assembly versions so we just need one of these
    # Doesn't actually matter which since the genome download will get the latest
    # version of the assembly automatically. Not all assemblies conform to the
    # ASMxxxxxx standard however...

    accessions = []
    assembly_re = re.compile(r'^(ASM[0-9]*)v([0-9])')
    for assembly in json_data:
        match = re.search(assembly_re, assembly['assembly_title'])
        if match:
            accessions.append(match.group(1))
            assembly['assembly'] = match.group(1)
            assembly['version'] = match.group(2)
        else:
            accessions.append(assembly['assembly_title'].split(maxsplit=1)[0])
            assembly['assembly'] = assembly['assembly_title'].split(maxsplit=1)[0]

    accessions=list(set(accessions))
    filtered_json_data = []

    for accession in accessions:
        for assembly in json_data:
            if assembly['assembly']==accession:
                filtered_json_data.append(assembly)
                break

    return filtered_json_data


def download_xml(local_file, uri):
    """
    Downloads xml file from remote site to specified directory

    Required parameters:
        local_file(pathlib::Path): Path for local copy of the xml
        uri(str): URI of xml document

    Returns:
        None
    """

    if not local_file.exists():

        result = make_request(uri)
        with open(local_file, "w", encoding='UTF-8') as fh:
            fh.writelines(result)

def get_xpath (tree, query):
    """
    Extracts desired field from XML tree with XPATH

    Required params:
        tree(etree): parsed xml tree
        query(str): Xpath query

    Returns:
        value(str): parsed value
    """
    tags = tree.xpath(query)

    if len(tags):
        return tags[0].text

    return None

def merge_fields(field1, field2):

    """ 
    Merges two fields if they are both present, with a ':' separator

    Required parameters: 
        field1(str): first field
        field2(str): second field

    Returns:
        merged(str): merged fields
    """

    merged = field2
    if field1 and field2:
        merged = ": ".join([field1, field2])

    return merged

def envo_lookup(term):
    """
    Looks up ENVO (environment ontology) terms to provide textual description

    These may be provided either as 'ENVO:xxxxxxxx' or just a numeric string where the prefix has not been provided
    
    Required params: 
        term(str): ENVO term
    
    Returns:
        text(str): label of ENVO term
    """
    if term:

        match = re.search(r'^(ENVO)?:?([0-9]+)$', term)
        if match:
            if (match.group(2) and not match.group(1)) or (match.group(1) and match.group(2)):
                term = f'ENVO_{match.group(2)}'
            elif not match.group(2):
                return None

            TERM_URI = f"{ENVO_URI}{term}"
            term_data = json.loads(make_request(TERM_URI))

            if 'label' in term_data.keys():
                return term_data['label']

    return term

def is_complete(local_file):
    """
    Determines whether an assembly is intact or not, based upon
    <= 5 scaffolds (to allow for plasmids)
    0 spanned gaps
    0 unspanned gaps

    Some assemblies are reporting 0 gaps when the are in hundreds of contigs, so 
    this is really an initial screen to save downloading clearly fragmented genomes
    but these will need checking once available locally

    Required params:
        local_file(pathlib.path) : Path to xml file

    Returns:
        state(bool): True if assembly is considered complete
    """

    state = False

    # set defaults which will not trigger detection if assembly
    # attributes missing from record
    spanned_gaps = 1
    unspanned_gaps = 1
    scaffold_count = False

    tree = etree.parse(local_file)

    # Some accessions are returning as not found so these need to be handled
    error = tree.xpath("/ErrorDetails/status")

    if error:
        print(f"{local_file} error status: {error[0].text}")
    else:

        #  Scaffold number and ghp info is stored within
        # ASSEMBLY_ATTRIBUTES as follow:

        #  <ASSEMBLY_ATTRIBUTE>
        #    <TAG>scaffold-count</TAG>
        #    <VALUE>1</VALUE>
        #  </ASSEMBLY_ATTRIBUTE>
        #
        #  <ASSEMBLY_ATTRIBUTE>
        #    <TAG>spanned-gaps</TAG>
        #    <VALUE>0</VALUE>
        #  </ASSEMBLY_ATTRIBUTE>
        #  <ASSEMBLY_ATTRIBUTE>
        #    <TAG>unspanned-gaps</TAG>
        #    <VALUE>0</VALUE>
        #  </ASSEMBLY_ATTRIBUTE>

        scaffold_count = get_xpath(tree, ".//ASSEMBLY_ATTRIBUTE/TAG[contains(text(),'scaffold-count')]/following-sibling::VALUE")
        spanned_gaps = get_xpath(tree, ".//ASSEMBLY_ATTRIBUTE/TAG[contains(text(),'spanned-gaps')]/following-sibling::VALUE")
        unspanned_gaps = get_xpath(tree, ".//ASSEMBLY_ATTRIBUTE/TAG[contains(text(),'unspanned-gaps')]/following-sibling::VALUE")

        if (
            scaffold_count
            and int(scaffold_count) <= 5
            and spanned_gaps == "0"
            and unspanned_gaps == "0"
        ):
            state = True

    return state


def query_biosample(biosample_acc, xml_dir, biosample_mapping):
    """
    Queries NCBI eutils for biosample metadata - EBI seems to be missing some entries which are present at the NCBI

    Required parameters:
        biosample_id(str): Biosample accession
        xml_dir(pathlib::Path): Path to biosample xmls
        biosample_mapping(dict): mapping of biosample id to NCBI internal ID


    Returns:
        metadata(dict): Parsed metadata values
        biosample_mapping(dict): biosample_mapping, possibly updated
    """

    fields = ["title", "strain", "culture_collection", "collection_date", "geographic location (region and locality)", "geo_loc_name",
               "lat_lon", "host", "isolation_source", "env_broad_scale", "env_local_scale"]

    metadata = dict.fromkeys(fields)

    biosample_id = biosample_mapping.get(biosample_acc)
    if not biosample_id:
        # eutils does not allow direct querying with the biosample accession - we first need to look up the ID
        uri = f"https://eutils.ncbi.nlm.nih.gov//entrez/eutils/esearch.fcgi?db=biosample&term={biosample_acc}&format=json"
        biosample_data = json.loads(make_request(uri))
        biosample_id = biosample_data["esearchresult"]["idlist"][0]
        biosample_mapping[biosample_acc] = biosample_id

    local_file = xml_dir / f"{biosample_id}.xml"
    # XML files are cached locally to ensure we can get a full set by rerunning due to vagaries of web services...
    if not local_file.exists():
        uri = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id={biosample_id}"
        download_xml(local_file, uri)

    with open(local_file, 'r', encoding='UTF-8') as fh:
        xml_doc = fh.read()

    tree = etree.fromstring(xml_doc.replace(' encoding="UTF-8"', ''))

    error = tree.xpath("/ErrorDetails/status")
    if error:
        print(f"{biosample_id} error status: {error[0].text}")
    else:
        query = "//BioSample/@submission_date"
        sub_date = tree.xpath(query)
        if sub_date:
            sub_date = datetime.strptime(tree.xpath(query)[0], '%Y-%m-%dT%H:%M:%S.%f')
        metadata['submission_date'] = sub_date

        for field in fields:
            # some field names may be in attribute_name or harmonized_name attributes
            query = f".//Attributes/Attribute[@attribute_name='{field}']"
            result = get_xpath(tree, query)

            if result:
                metadata[field] = result
            else:
                query = f".//Attributes/Attribute[@harmonized_name='{field}']"
                metadata[field] = get_xpath(tree, query)

    return metadata, biosample_mapping

def extract_metadata(genome, local_xml, xml_dir, biosample_mapping):
    """
    Extracts some metadata from the xml file and also via query_biosample
    
    Required parameters: 
        genome: (pathlib:Path): Path to genome record
        local_xml(pathlib::Path): Path to local xml file for assembly
        xml_dir(pathlib::Path): Path to store biosample xml file
        biosample_mapping(dict): Mapping of biosample_ids to NCBI internal IDs
    
    Returns:
        metadata(dict) - dictionary of collected values
        biosample_mapping(dict) - mapping of biosample id to internal NCBI id
    """

    tree = etree.parse(local_xml)

    metadata = dict.fromkeys(['study_id', 'biosample_id'])

    # Some accessions are returning as not found so these need to be handled
    error = tree.xpath("/ErrorDetails/status")
    if error:
        print(f"{local_xml} error status: {error[0].text}")
    else:
        tags = tree.xpath( ".//STUDY_REF/IDENTIFIERS/PRIMARY_ID")
        if len(tags):
            metadata['study_id'] = tags[0].text

        tags = tree.xpath( ".//SAMPLE_REF/IDENTIFIERS/PRIMARY_ID")
        if len(tags):
            metadata['biosample_id'] = tags[0].text

    biosample_metadata, biosample_mapping = query_biosample(metadata['biosample_id'], xml_dir, biosample_mapping)
    metadata = biosample_metadata | metadata

    # add additional metadate from genome record
    metadata = extract_embl_metadata(metadata, genome)

    return metadata, biosample_mapping

def extract_embl_metadata(metadata, genome):
    """
    Attempt to populate metadata fields not present in biosample metadata from embl record

    Required parameters:
        metadata(dict): metadata dictionary
        genome(pathlib::Path): Path to downloaded genome record
    
    Returns:
        metadata(dict): metadata dictionary
    """

    if genome.stat().st_size > 20: # size of an empty gzipped record
        with gzip.open(genome, "rt") as embl_fh:
            record = list(SeqIO.parse(embl_fh, format="embl"))[0]

            if not metadata['title']:
                metadata["title"] = clean_description(record.description)

            if not metadata["isolation_source"] or not metadata["strain"]:
                for feature in record.features:
                    if feature.type == "source":
                        source = feature.qualifiers.get("isolation_source")
                        if source is not None:
                            metadata["isolation_source"] = source[0]

                        isolate = feature.qualifiers.get("isolate")
                        if isolate is not None and not metadata["strain"] :
                            metadata["strain"] = isolate[0]

                        break

    return metadata

def download_assembly(accession, local_genome):

    """
    Retrieves genome assembly from ENA

    Required parameters:
        accession(str): Accession to retrieve
        local_genome(pathlib.Path): Path to store file locally

    Returns:
        ok(bool): True on successful download
    """

    ok = True

    if not local_genome.exists():

        uri = f"{ENA_URI}browser/api/embl/{accession}?download=true&gzip=true"

        response = requests.get(uri, timeout=30)
        if response.status_code == 200:
            with open(local_genome, "wb") as fh:
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        fh.write(chunk)
        else:
            print(f"Error retrieving {accession}")
            print(response.status_code)
            ok = False

    return ok


def fasta_convert(accession, genome_dir, fasta_dir):
    """
    Converts an embl format records to fasta

    Required parameters:
        accession(str): Accession of entry to convert
        genome_dir(Path): Path to genomes directory
        fasta_dir(Path): Path to fasta directory

    Returns:
        scaffolds(int): number of scaffolds in assembly
        seq_ids(list): list of sequence ids in assembly
    """

    records = []
    seq_ids = []

    genome = genome_dir / f"{accession}.embl.gz"
    fasta_dir = fasta_dir / "genomes"
    fasta_dir.mkdir(exist_ok=True)

    fasta_path = fasta_dir / f"{accession}.fasta"

    with gzip.open(genome, "rt") as embl_fh:
        for record in SeqIO.parse(embl_fh, format="embl"):
            seq_ids.append(record.id)

            record.id = f"lcl|{accession}|{record.id}"
            records.append(record)

    if records and not fasta_path.exists():  # for empty embl records from bad accessions
        with open(fasta_path, "w", encoding='UTF-8') as fasta_fh:
            SeqIO.write(records, fasta_fh, format="fasta")

    return len(records), seq_ids


def main():
    """Main process"""


    parser = argparse.ArgumentParser(
        prog = __file__, 
        description="Downloads, fasta converts and optionally blast indexes genome records for organism"
    )
    parser.add_argument('-b', '--blast', action='store_true', required=False, help="Generate blast index")
    args = parser.parse_args()

    species_name = get_taxa_name(NCBI_TAXID)
    metadata_list = []

    # Locations for storing various data types
    ena_metadata_dir = Path(__file__).parents[1] / Path("data/full/metadata/ena")
    biosample_metadata_dir = Path(__file__).parents[1] / Path("data/full/metadata/biosample")
    genome_dir = Path(__file__).parents[1] / Path("data/full/genomes")
    fasta_dir = Path(__file__).parents[1] / Path("data/full/fasta")
    blast_dir = Path(__file__).parents[1] / Path("data/full/blast_db")
    output_dir = Path(__file__).parents[1] / Path("data/full/outputs")
    biosample_mapping_file = biosample_metadata_dir / 'biosamples.json'
    id_mapping_file = Path(__file__).parents[1] / Path("data/full/id_mapping.json")

    for data_dir in (ena_metadata_dir, biosample_metadata_dir, genome_dir, fasta_dir, blast_dir, output_dir):
        data_dir.mkdir(exist_ok=True, parents=True)

    # Obtain total list of assemblies available
    genome_info = search_available()

    # Mapping of seq ids -> GCA number/strain name
    id_mapping = {}

    # Read biosample ID mapping if it exists, otherwise create empty dict for it....
    biosample_mapping = {}

    if biosample_mapping_file.exists():
        with open(biosample_mapping_file, 'r', encoding='UTF-8') as fh:
            biosample_mapping = json.load(fh)

    for assembly in tqdm(genome_info, desc="Download", leave=None):

        accession = assembly.get("accession")
        local_ena_xml = ena_metadata_dir / Path(f"{accession}.xml")
        local_genome = genome_dir / Path(f"{accession}.embl.gz")

        # Download SRA assembly XML...
        download_xml(local_ena_xml, f"{ENA_URI}browser/api/xml/{accession}")
        # Check for genome completeness
        complete = is_complete(local_ena_xml)

        if complete:
            ok = download_assembly(accession, local_genome)
            if ok:
                # Collect some metadata
                metadata, biosample_mapping = extract_metadata(local_genome, local_ena_xml, biosample_metadata_dir, biosample_mapping)
                metadata['accession'] = accession
                scaffolds, seq_ids = fasta_convert(accession, genome_dir, fasta_dir)
                metadata['scaffolds'] = scaffolds
                metadata_list.append(metadata)
                for seq_id in seq_ids:
                    id_mapping[seq_id] = {
                        'accession': accession,
                        'strain': metadata['strain']
                    }
            else:
                print(f"{accession} could not be successfully downloaded")

    #if args.blast:
    #    blast_index(fasta_dir, blast_dir, species_name, NCBI_TAXID, 'nucl')

    with open(biosample_mapping_file, 'w', encoding='UTF-8') as out_fh:
        json.dump(biosample_mapping, out_fh)

    with open(id_mapping_file, 'w', encoding='UTF-8') as out_fh:
        json.dump(id_mapping, out_fh)

    summary = pd.DataFrame(metadata_list)
    summary = summary[summary.scaffolds != 0]

    # remove the various different ways people provide metadata without providing metadata...
    for phrase in ["[Nn]ot [Aa]pplicable", "NA", "missing", "[Uu]nknown", "[Nn]ot collected"]:
        summary.replace(phrase, "", regex=True, inplace=True)

    # Look up any ENVO terms which are used to provide environmental metadata
    for field in ["isolation_source", "env_broad_scale", "env_local_scale"]:
        summary[field] = summary.apply(lambda x: envo_lookup(x[field]), axis=1)

    summary['location'] = summary.apply(lambda x: merge_fields(x["geographic location (region and locality)"], x["geo_loc_name"]), axis=1)
    summary['env_context'] = summary.apply(lambda x: merge_fields(x["env_broad_scale"], x["env_local_scale"]), axis=1)
    summary.drop(["geographic location (region and locality)", "geo_loc_name", "env_broad_scale", "env_local_scale"], axis=1, inplace=True)

    cols = ["accession", "study_id", "biosample_id", "title", "strain", "culture_collection",
            "host", "isolation_source", "env_context", "location", "collection_date", "scaffolds"]
    summary = summary[cols]

    summary.rename({'accession': 'Accession', 'study_id': 'Study ID', 'biosample_id': 'Biosample ID', 'title': 'Title',
                    'strain': 'Strain', 'culture_collection': 'Culture collection', 'host': 'Host', 
                    'isolation_source': 'Isolation source', 'env_context': 'Environmental context', 
                    'location': 'Location', 'collection_date': 'Collection date', 'scaffolds': 'Scaffolds' }, axis=1, inplace=True)

    date = datetime.today().strftime("%d-%m-%Y")
    outfile = f"data/full/{species_name.replace(' ', '_')}_complete_genomes_{date}.txt"
    summary.to_csv(outfile, sep="\t", index=False, header=True)

if __name__ == "__main__":
    main()
