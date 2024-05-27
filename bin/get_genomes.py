#!/bin/env python

import gzip
import requests

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from json import loads
from lxml import etree
from pathlib import Path
from tqdm import tqdm

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

    uri = f"{ENA_URI}portal/api/search?result=assembly&query=tax_tree(1423)&fields=assembly_title,tax_id&format=json"
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

    tree = etree.parse(local_file)

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
    spanned_gaps = tags[0].text

    tags = tree.xpath(".//ASSEMBLY_ATTRIBUTE/TAG[contains(text(),'unspanned-gaps')]/following-sibling::VALUE")
    unspanned_gaps = tags[0].text

    state = False

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

    required_genes = ['comP']
    gene_seqs = dict()
    accession = str(local_genome.stem).replace('.embl','')

    with gzip.open(local_genome, 'rt') as fh:
        for record in SeqIO.parse(fh, format = 'embl'):
            for feature in record.features:

                # Handle genes by extracting subsequence from record
                if feature.type == 'gene':
                    if 'gene' in feature.qualifiers.keys() and feature.qualifiers.get('gene')[0] in required_genes:

                        gene = feature.extract(record)
                        id = f'{accession}_{feature.qualifiers.get('gene')[0]}_gene'
                        gene.id = id
                        gene_seqs[id] = gene

                # CDS feature contain the protein sequence as a 'translation' qualifier
                elif feature.type == 'CDS':

                    if 'gene' in feature.qualifiers.keys() and feature.qualifiers.get('gene')[0] in required_genes:
                        if 'translation' in feature.qualifiers.keys():

                            prot_seq = feature.qualifiers.get('translation')[0]
                            id = f'{accession}_{feature.qualifiers.get('gene')[0]}_protein'
                            protein = SeqRecord(Seq(prot_seq), id = id)

                            gene_seqs[id] = protein

    return(gene_seqs)


def main():

    xml_dir = Path(__file__).parents[1] / Path('xml')
    genome_dir = Path(__file__).parents[1] / Path('genomes')

    xml_dir.mkdir(exist_ok = True)
    genome_dir.mkdir(exist_ok = True)

    genome_info = search_available()

    all_seqs = list()

    for assembly in tqdm(genome_info):

        accession = assembly.get('accession')

        local_xml = xml_dir / Path(f"{accession}.xml")
        local_genome = genome_dir / Path(f"{accession}.embl.gz")

        download_assembly_xml(local_xml, accession)
        complete = is_complete(local_xml)

        if complete:
            download_assembly(accession, local_genome)
            gene_seqs = get_gene_sequences(local_genome)
            all_seqs.append(gene_seqs)

        break

if __name__ == "__main__":
    main()