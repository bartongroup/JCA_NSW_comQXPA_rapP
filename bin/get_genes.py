#!/bin/env python

import gzip
import pandas as pd
import re
import requests

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import ChainMap
from json import loads
from lxml import etree
from pathlib import Path
import subprocess
from tqdm import tqdm

from pprint import pprint

NCBI_TAXID     = '1423'
REQUIRED_GENES = ['comP', 'walK', 'yycG', 'desK', 'yocF']
ENA_URI="https://www.ebi.ac.uk/ena/"

"""
Identifies all complete B.subtilis genomes available via ENA and downloads them
"""

def make_request(uri):

    """
    Makes HTTP request and returns result body

    Required params:
        uri(str): URI for request

    Returns:
        text of response
    """
    try:
        r = requests.get(uri)
    except Exception as e:
        print(f'Error requesting {uri}: {e}')
    
    return(r.text)

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
    json = loads(result)

    return(json)

def download_assembly_xml(local_file, accession):

    """
    Retrieves details on an assembly from EBI

    Required parameters:
        local_file(pathlib::Path): Path for local copy of the xml
        accession(str): assembly accession

    Returns:
        None
    """

    if not local_file.exists():

        uri = f"{ENA_URI}browser/api/xml/{accession}"
        result = make_request(uri)
        with open(local_file, 'w') as fh:
            fh.writelines(result)

def is_complete(local_file):

    """
    Determines whether an assembly is intact or not, based upon
    <= 5 scaffolds (to allow for plasmids)
    0 spanned gaps
    0 unspanned gaps

    Some assemblies are reporting 0 gaps when the are in hundreds of contigs...
    
    Required params: 
        local_file(pathlib.path) : Path to xml file

    Returns:
        state(bool): True if assembly is considered complete 
    """

    state = False

    # set defaults which will not trigger detection if assembly attributes missing from record
    spanned_gaps = 1
    unspanned_gaps = 1
    scaffold_count = False

    tree = etree.parse(local_file)

    # Some accessions are returning as not found so these need to be handled
    error = tree.xpath("/ErrorDetails/status")

    if error:
        print(f"{local_file} error status: {error[0].text}")
    else:
        #  Scaffold number and ghp info is stored within ASSEMBLY_ATTRIBUTES as follow:
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

        tags = tree.xpath(".//ASSEMBLY_ATTRIBUTE/TAG[contains(text(),'scaffold-count')]/following-sibling::VALUE")
        if len(tags):
            scaffold_count = tags[0].text

        tags = tree.xpath(".//ASSEMBLY_ATTRIBUTE/TAG[contains(text(),'spanned-gaps')]/following-sibling::VALUE")
        if len(tags):
            spanned_gaps = tags[0].text

        tags = tree.xpath(".//ASSEMBLY_ATTRIBUTE/TAG[contains(text(),'unspanned-gaps')]/following-sibling::VALUE")
        if len(tags):
            unspanned_gaps = tags[0].text

        if scaffold_count and int(scaffold_count) <= 5 and spanned_gaps == '0' and unspanned_gaps == '0':
            state = True

    return(state)

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

        response = requests.get(uri)
        if response.status_code == "200":
            with open(local_genome, 'wb') as fh:
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        fh.write(chunk)
        else:
            print(f"Error retrieving {accession}")
            ok = False

    return(ok)

def clean_description(description):
    """
    tidies up variously formatted descriptions into something reasonably consistent

    Required params:
        description(str): record description
    
    Returns:
        description(str): cleaned description
    """

    description = re.sub(r'[, ]*complete genome[.]*', '', description)
    description = re.sub(r'[, ]*genome assembly[,:_ A-Za-z0-9]*', '', description)
    description = re.sub(r'chromosome[ ,]*', '', description)
    description = re.sub('whole genome shotgun sequence.','', description)
    description = re.sub(r'NODE[A-Za-z0-9_,]*', '', description)
    description = re.sub(r'[ ,]*$','', description)

    return(description)


def fasta_convert(genome_dir, fasta_dir):
    """
    Converts all embl format records to fasta 

    Required parameters:
        genome_dir(Path): Path to genomes directory
        fasta_dir(Path): Path to fasta directory

    Returns:
        summary(pd.DataFrame): Data frame containing metadata on assemblies
    """

    genomes = genome_dir.glob("*.gz")
    metadata_list = list()

    for genome in tqdm(genomes, desc = 'Reformat', leave = None):

        metadata = dict()
        records  = list()

        accession  = genome.name.replace(".embl.gz","")
        fasta_file = f"{accession}.fasta"
        fasta_path = fasta_dir / Path(fasta_file)

        with gzip.open(genome, 'rt') as embl_fh:
            for record in SeqIO.parse(embl_fh, format = 'embl'):

                if not metadata.get('accession'):

                    metadata['accession'] = accession
                    metadata['isolate'] = clean_description(record.description)

                    for feature in record.features:
                        if feature.type == 'source':

                            source = feature.qualifiers.get('isolation_source')
                            if source is not None:
                                metadata['source'] = source[0]
                record.id = f"lcl|{record.id}"
                records.append(record)

        with open(fasta_path, 'w') as fasta_fh:
            SeqIO.write(records, fasta_fh, format = 'fasta')

        metadata['scaffolds'] = len(records)
        metadata_list.append(metadata)

    summary = pd.DataFrame(metadata_list)
    summary.to_excel('summary.xlsx', sheet_name = 'Complete genomes', index = False, header = True)

    return(summary)

def blast_index(fasta_dir):

    """
    Creates a blast index for each fasta file contained within directory

    Required parameters:
        fasta_dir(Path): Location of fasta files

    Returns:
        None
    """

    fasta_files = fasta_dir.glob('*.fasta')
    for fasta_file in tqdm(fasta_files, desc = 'Blast index'):

        accession = fasta_file.name.replace('.fasta','')

        cmd = ["makeblastdb",  "-in",  fasta_file, "-dbtype", "nucl", "-taxid", NCBI_TAXID, "-title", accession]
        try:
            result = subprocess.run(cmd, check = True, capture_output = True)
        except subprocess.CalledProcessError as e:
            print(f"Index failed: {accession}")
            print(result.stdout)
            print(result.stderr)

def get_gene_sequences(local_genome):

    """
    Parses the downloaded genome record to extract the required sequences

    required parameters:
        local_genome(libpath.Path): Path to local copy of record

    returns:
        gene_seqs(dict): dictionary containing gene and amino acid sequences 
                         as BioSeq objects, keyed on seq ID
    """

    gene_seqs = dict()
    accession = str(local_genome.stem).replace('.embl','')

    with gzip.open(local_genome, 'rt') as fh:
        for record in SeqIO.parse(fh, format = 'embl'):
            organism = record.annotations['organism']

            for feature in record.features:

                # Handle genes by extracting subsequence from record
                if feature.type == 'gene':
                    if 'gene' in feature.qualifiers.keys() and feature.qualifiers.get('gene')[0] in REQUIRED_GENES:

                        gene = feature.extract(record)
                        gene_id = feature.qualifiers.get('gene')[0]

                        id = f'{accession}_{gene_id}_gene'

                        gene.id = id
                        gene.description = f"{organism} {gene_id} gene"

                        gene_seqs[id] = gene

                # CDS feature contain the protein sequence as a 'translation' qualifier
                elif feature.type == 'CDS':

                    if 'gene' in feature.qualifiers.keys() and feature.qualifiers.get('gene')[0] in REQUIRED_GENES:

                        if 'translation' in feature.qualifiers.keys():

                            gene_id = feature.qualifiers.get('gene')[0]
                            prot_seq = feature.qualifiers.get('translation')[0]
                            id = f'{accession}_{gene_id}_protein'

                            protein = SeqRecord(Seq(prot_seq), 
                                                id = id, 
                                                description = f"{organism} {gene_id} protein")

                            gene_seqs[id] = protein

    return(gene_seqs)

def write_outputs(all_seqs, output_dir):

    """
    Writes fasta formatted records for each gene/protein

    Required params:
        all_seqs(dict): dicts combing results returned from get_gene_sequences()
        output_dir(pathlib.Path): Location to write results

    Returns:
        None
    """

    separated_seqs = dict()

    for gene in REQUIRED_GENES:
        for seq_id in all_seqs.keys():
            record = all_seqs[seq_id]

            if gene in seq_id:

                match = re.search(r'_(gene|protein)$', seq_id)
                if match:
                    seq_type = match.group(1)
                else:
                    raise Exception (f'Failed to match sequence type: {seq_id}')

                record.id = seq_id.replace(f"_{gene}_{seq_type}", "")
                if f"{gene}_{seq_type}" in separated_seqs:
                    separated_seqs[f"{gene}_{seq_type}"].append(record)
                else:
                    separated_seqs[f"{gene}_{seq_type}"] = [record]

    for key in separated_seqs.keys():
        gene,seq_type = key.split('_')
        with open(f'outputs/{gene}.{seq_type}.fa', 'w') as fh:
            SeqIO.write(separated_seqs[key], fh, 'fasta')

def main():

    # Locations for storing various data types
    xml_dir    = Path(__file__).parents[1] / Path('xml')
    genome_dir = Path(__file__).parents[1] / Path('genomes')
    fasta_dir  = Path(__file__).parents[1] / Path('fasta')
    output_dir = Path(__file__).parents[1] / Path('outputs')

    xml_dir.mkdir(exist_ok = True)
    genome_dir.mkdir(exist_ok = True)
    fasta_dir.mkdir(exist_ok = True)
    output_dir.mkdir(exist_ok = True)

    # list for storing combined retrieved sequences
    all_seqs = list()

    # Obtain total list of assemblies available
    genome_info = search_available()

    for assembly in tqdm(genome_info, desc = 'Download', leave = None):

        accession = assembly.get('accession')

        local_xml    = xml_dir    / Path(f"{accession}.xml")
        local_genome = genome_dir / Path(f"{accession}.embl.gz")

        # Download SRA assembly XML...
        download_assembly_xml(local_xml, accession)
        # Check for completeness
        complete = is_complete(local_xml)

        if complete:
            # Download EMBL formatted record
            ok = download_assembly(accession, local_genome)
            if ok:
                # Parse required gene/protein sequences
                gene_seqs = get_gene_sequences(local_genome)
                # Add results to overall list
                all_seqs.append(gene_seqs)

    # Convert assemblies to fasta format, which also provides
    # a summary pandas dataframe of the data
    summary = fasta_convert(genome_dir, fasta_dir)
    blast_index(fasta_dir)
            
    # Merge dicts into one giant dict...
    all_seqs = dict(ChainMap(*all_seqs))

    #...and create the outputs
    write_outputs(all_seqs, output_dir)

if __name__ == "__main__":
    main()