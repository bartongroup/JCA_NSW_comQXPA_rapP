#!/bin/env python

"""
Creates a database of all 'complete' genome sequences for a TAXID, reformats
them to fasta, while capturing metadata from biosample and the genome record
"""

import gzip
from pathlib import Path
from datetime import datetime
import json

from Bio import SeqIO
from lxml import etree
from tqdm import tqdm
import pandas as pd
import requests

from common import (
    get_taxa_name,
    clean_description,
    make_request,
)

NCBI_TAXID = "1423"
ENA_URI = "https://www.ebi.ac.uk/ena/"


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

    return json_data


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
    else:
        return None


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


def query_biosample(biosample_acc, xml_dir):
    """
    Queries NCBI eutils for biosample metadata - EBI seems to be missing some entries which are present at the NCBI

    Required parameters:
        biosample_id(str): Biosample accession

    Returns:
        metadata(dict): Parsed metadata values
    """

    fields = ["strain", "culture_collection", "collection_date", "geo_loc_name", "host", "isolation_source"]

    metadata = dict.fromkeys(fields)

    # eutils does not allow direct querying with the biosample accession - we first need to look up the ID
    uri = f"https://eutils.ncbi.nlm.nih.gov//entrez/eutils/esearch.fcgi?db=biosample&term={biosample_acc}&format=json"
    biosample_data = json.loads(make_request(uri))
    biosample_id = biosample_data["esearchresult"]["idlist"][0]

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
        query = ".//Biosample/@submission_date"
        metadata['submission_date'] = get_xpath(tree, query)

        for field in fields:
            # some field names may be in attribute_name or harmonized_name attributes
            query = f".//Attributes/Attribute[@attribute_name='{field}']"
            result = get_xpath(tree, query)

            if result:
                metadata[field] = result
            else:
                query = f".//Attributes/Attribute[@harmonized_name='{field}']"
                metadata[field] = get_xpath(tree, query)

        if not metadata['geo_loc_name']:
            # we have some misformed tags...in our own submissions...oops...
            query = ".//Attributes/Attribute[@attribute_name='geographic location (region and locality']"
            metadata['geo_loc_name'] = get_xpath(tree, query)

    return metadata

def extract_metadata(genome, local_xml, xml_dir):
    """
    Extracts some metadata from the xml file and also via query_biosample
    
    Required parameters: 
        local_xml(pathlib::Path) - Path to local xml file for assembly
        xml_dir(pathlib::Path) - Path to store biosample xml file
    
    Returns:
        metadata(dict) - dictionary of collected values
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

    biosample_metadata = query_biosample(metadata['biosample_id'], xml_dir)
    metadata = metadata | biosample_metadata

    # add additional metadate from genome record
    metadata = extract_embl_metadata(metadata, genome)

    return metadata

def extract_embl_metadata(metadata, genome):
    """
    Attempt to populate metadata fields not present in biosample metadata from embl record

    Required parameters:
        metadata(dict): metadata dictionary
        genome(pathlib::Path): Path to downloaded genome record
    
    Returns:
        metadata(dict): metadata dictionary
    """

    # 
    if genome.stat().st_size > 20: # size of an empty gzipped record
        with gzip.open(genome, "rt") as embl_fh:
            record = list(SeqIO.parse(embl_fh, format="embl"))[0]

            metadata["title"] = clean_description(record.description)

            if not metadata["isolation_source"]:
                for feature in record.features:
                    if feature.type == "source":
                        source = feature.qualifiers.get("isolation_source")
                        if source is not None:
                            metadata["isolation_source"] = source[0]
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
        None
    """

    records = []
    genome = genome_dir / f"{accession}.embl.gz"

    fasta_path = fasta_dir / f"{accession}.fasta"

    with gzip.open(genome, "rt") as embl_fh:
        for record in SeqIO.parse(embl_fh, format="embl"):

            record.id = f"lcl|{record.id}"
            records.append(record)

    if records:  # for empty embl records from bad accessions
        with open(fasta_path, "w", encoding='UTF-8') as fasta_fh:
            SeqIO.write(records, fasta_fh, format="fasta")

    return len(records)


def main():
    """Main process"""
    species_name = get_taxa_name(NCBI_TAXID)
    metadata_list = []

    # Locations for storing various data types
    ena_xml_dir = Path(__file__).parents[1] / Path("data/full/xml/ena")
    biosample_xml_dir = Path(__file__).parents[1] / Path("data/full/xml/biosample")
    genome_dir = Path(__file__).parents[1] / Path("data/full/genomes")
    fasta_dir = Path(__file__).parents[1] / Path("data/full/fasta")
    blast_dir = Path(__file__).parents[1] / Path("data/full/blast_db")
    output_dir = Path(__file__).parents[1] / Path("data/full/outputs")

    for data_dir in (ena_xml_dir, biosample_xml_dir, genome_dir, fasta_dir, blast_dir, output_dir):
        data_dir.mkdir(exist_ok=True, parents=True)

    # Obtain total list of assemblies available
    genome_info = search_available()

    for assembly in tqdm(genome_info, desc="Download", leave=None):

        accession = assembly.get("accession")

        local_ena_xml = ena_xml_dir / Path(f"{accession}.xml")
        local_genome = genome_dir / Path(f"{accession}.embl.gz")

        # Download SRA assembly XML...
        download_xml(local_ena_xml, f"{ENA_URI}browser/api/xml/{accession}")
        # Check for genome completeness
        complete = is_complete(local_ena_xml)

        if complete:
            ok = download_assembly(accession, local_genome)
            if ok:
                # Collect some metadata
                metadata = extract_metadata(local_genome, local_ena_xml, biosample_xml_dir)
                metadata['accession'] = accession

                scaffolds = fasta_convert(accession, genome_dir, fasta_dir)
                metadata['scaffolds'] = scaffolds
                metadata_list.append(metadata)
            else:
                print(f"{accession} could not be successfully downloaded")

    summary = pd.DataFrame(metadata_list)
    date = datetime.today().strftime("%d-%m-%Y")
    outfile = f"{species_name.replace(' ', '_')}_complete_genomes_{date}.txt"
    summary = summary[summary.scaffolds != 0]

    summary.to_csv(outfile, sep="\t", index=False, header=True)

    summary = summary[['accession', 'isolate']]
    summary.to_csv("strains.txt", sep="\t", header=True, index=False)


if __name__ == "__main__":
    main()
