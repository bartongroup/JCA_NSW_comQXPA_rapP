#!/bin/env python

import gzip
import pandas as pd
import re
import requests

from Bio import SeqIO
from datetime import datetime
from json import loads
from lxml import etree
from tqdm import tqdm
from pathlib import Path

from pprint import pprint

from common import blast_index, get_taxa_name, clean_description, clean_strain, make_request

NCBI_TAXID = '1423'
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
        if response.status_code == 200:
            with open(local_genome, 'wb') as fh:
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        fh.write(chunk)
        else:
            print(f"Error retrieving {accession}")
            print(response.status_code)
            ok = False

    return(ok)


def fasta_convert(genome_dir, fasta_dir, species_name):

    """
    Converts all embl format records to fasta 

    Required parameters:
        genome_dir(Path): Path to genomes directory
        fasta_dir(Path): Path to fasta directory
        species_name(str): Name of species

    Returns:
        None
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
                    metadata['clean_strain'] = clean_strain(metadata['isolate'])

                    for feature in record.features:
                        if feature.type == 'source':

                            source = feature.qualifiers.get('isolation_source')
                            if source is not None:
                                metadata['source'] = source[0]
                            
                record.id = f"lcl|{record.id}"
                records.append(record)

        if len(records): # for empty embl records
            with open(fasta_path, 'w') as fasta_fh:
                SeqIO.write(records, fasta_fh, format = 'fasta')

            metadata['scaffolds'] = len(records)
            metadata_list.append(metadata)

    summary = pd.DataFrame(metadata_list)
    date = datetime.today().strftime('%d-%m-%Y')
    outfile = f"{species_name.replace(' ','_')}_complete_genomes_{date}.xlsx"
    summary = summary[summary.scaffolds != 0]
    summary.to_excel(outfile, sheet_name = 'Complete genomes', index = False, header = True)
    summary.drop(['source', 'scaffolds', 'clean_strain'], inplace = True, axis = 1)
    summary.to_csv('strains.txt', sep="\t", header = True, index = False)

    return(summary)

def main():

    species_name = get_taxa_name(NCBI_TAXID)

    # Locations for storing various data types
    xml_dir    = Path(__file__).parents[1] / Path('xml')
    genome_dir = Path(__file__).parents[1] / Path('genomes')
    fasta_dir  = Path(__file__).parents[1] / Path('fasta')
    blast_dir  = Path(__file__).parents[1] / Path('blast_db')
    output_dir = Path(__file__).parents[1] / Path('outputs')

    xml_dir.mkdir(exist_ok = True)
    genome_dir.mkdir(exist_ok = True)
    fasta_dir.mkdir(exist_ok = True)
    blast_dir.mkdir(exist_ok = True)
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
            #if ok:
                # Parse required gene/protein sequences
            #    gene_seqs = get_gene_sequences(local_genome)
            #    # Add results to overall list
            #    all_seqs.append(gene_seqs)

    # Convert assemblies to fasta format, which also provides
    # a summary pandas dataframe of the data
    fasta_convert(genome_dir, fasta_dir, species_name)
    blast_index(fasta_dir, blast_dir, species_name)

if __name__ == "__main__":
    main()
