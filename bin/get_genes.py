#!/bin/env python

import gzip
import re
import requests

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import ChainMap
from json import loads
from lxml import etree
from pathlib import Path
from tqdm import tqdm

from pprint import pprint

NCBI_TAXID     = '1423'
REQUIRED_GENES = ['comP']

"""
Identifies all complete B.subtilis genomes available via ENA and downloads them
"""

ENA_URI="https://www.ebi.ac.uk/ena/"

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
    Determines whether an assembly is intact or not
    
    Required params: 
        local_file(pathlib.path) : Path to xml file

    Returns:
        state(bool): True if assembly contains no gaps
    """

    state = False

    # set defaults which will not trigger detection if assembly attributes missing from record
    spanned_gaps = 1
    unspanned_gaps = 1

    tree = etree.parse(local_file)

    # Some accessions are returning as not found so these need to be handled
    error = tree.xpath("/ErrorDetails/status")

    if error:
        print(f"{local_file} error status: {error[0].text}")
    else:
        #  Gap info is stored within ASSEMBLY_ATTRIBUTES as follow:
        #
        #  <ASSEMBLY_ATTRIBUTE>
        #   <TAG>spanned-gaps</TAG>
        #   <VALUE>0</VALUE>
        #  </ASSEMBLY_ATTRIBUTE>
        #  <ASSEMBLY_ATTRIBUTE>
        #    <TAG>unspanned-gaps</TAG>
        #    <VALUE>0</VALUE>
        #  </ASSEMBLY_ATTRIBUTE>

        tags = tree.xpath(".//ASSEMBLY_ATTRIBUTE/TAG[contains(text(),'spanned-gaps')]/following-sibling::VALUE")
        if len(tags):
            spanned_gaps = tags[0].text

        tags = tree.xpath(".//ASSEMBLY_ATTRIBUTE/TAG[contains(text(),'unspanned-gaps')]/following-sibling::VALUE")
        if len(tags):
            unspanned_gaps = tags[0].text


        if spanned_gaps == '0' and unspanned_gaps == '0':
            state = True

    return(state)

def download_assembly(accession, local_genome):

    """
    Retrieves genome assembly from ENA

    Required parameters:
        accession(str): Accession to retrieve
        local_genome(pathlib.Path): Path to store file locally
    """

    if not local_genome.exists():

        uri = f"{ENA_URI}browser/api/embl/{accession}?download=true&gzip=true"

        try:
            response = requests.get(uri)
        except Exception as e:
            print(f'Error requesting {uri}: {e}')

        with open(local_genome, 'wb') as fh:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    fh.write(chunk)

def get_gene_sequences(local_genome):

    """
    Parses the downloaded genome record to extract the comP and xxxxx[TODO: TBD] sequences

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
    output_dir = Path(__file__).parents[1] / Path('outputs')

    xml_dir.mkdir(exist_ok = True)
    genome_dir.mkdir(exist_ok = True)
    output_dir.mkdir(exist_ok = True)

    # list for storing combined retrieved sequences
    all_seqs = list()

    # Obtain total list of assemblies available
    genome_info = search_available()
    complete_count = 0

    for assembly in tqdm(genome_info):

        accession = assembly.get('accession')

        local_xml    = xml_dir    / Path(f"{accession}.xml")
        local_genome = genome_dir / Path(f"{accession}.embl.gz")

        # Download SRA assembly XML...
        download_assembly_xml(local_xml, accession)
        # Check for completeness
        complete = is_complete(local_xml)

        if complete:
            complete_count += 1
            # Download EMBL formatted record
            download_assembly(accession, local_genome)
            # Parse required gene/protein sequences
            gene_seqs = get_gene_sequences(local_genome)
            # Add results to overall list
            all_seqs.append(gene_seqs)
            
    # Merge dicts into one giant dict...
    all_seqs = dict(ChainMap(*all_seqs))

    #...and create the outputs
    write_outputs(all_seqs, output_dir)

    print(f'Complete genomes: {complete_count}')

if __name__ == "__main__":
    main()