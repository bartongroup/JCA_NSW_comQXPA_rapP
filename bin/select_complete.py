#!/bin/env python 

''' 
Selects genomes based upon BUSCO completion criteria
Looking for 98% completeness or better...
'''
THRESHOLD  = 98
NCBI_TAXID = '1423'

from json import load
from pathlib import Path
from shutil import copy
from tqdm import tqdm
from common import blast_index, get_taxa_name

def check_completeness():

    """
    Checks completeness of each genome based on parsing 'Complete percentage'
    score from BUSCO short summary json file. The threshold for comparison is
    defined in the THRESHOLD constant

    Required params:
        None

    Returns:
        retain(list): Accessions passing threshold
    """

    results = list(Path('busco').glob('*/short_summary.specific*json'))
    retain = list()

    for result in tqdm(results):

        file = result.name
        accession = file.replace('short_summary.specific.bacillales_odb10.', '')
        accession = accession.replace('.json','')

        with open(result) as fh:
            dat = load(fh)
            complete = dat.get('results').get('Complete percentage')

            if float(complete) > THRESHOLD:
                retain.append(accession)

    print(f"{len(retain)}/{len(results)} genomes exceed completeness threshold of {THRESHOLD}%")

    return(retain)

def create_good_set(species_name, retain):

    """
    Creates a subset of the databases flagged for retention

    Required params:
        species_name(str): Name of species
        retain(list): Accessions of sequences to retain

    Returns:
        none
    """

    annotated_genomes = Path('annotations')
    original_fastas = Path('fasta')
    complete_genomes = Path('complete/genomes')
    complete_fasta = Path('complete/fasta')
    complete_blast = Path('complete/blast')

    complete_genomes.mkdir(exist_ok = True, parents = True)
    complete_fasta.mkdir(exist_ok = True, parents = True)
    complete_blast.mkdir(exist_ok = True, parents = True)

    for accession in tqdm(retain):

        source_genome = annotated_genomes / accession / f'{accession}.embl'
        source_fasta =original_fastas / f'{accession}.fasta'
        destination_genome = complete_genomes / f'{accession}.embl'
        destination_fasta = complete_fasta / f'{accession}.fasta'
        
        copy(source_genome, destination_genome)
        copy(source_fasta, destination_fasta)
    
    blast_index(complete_fasta, complete_blast, species = species_name)


def main():

    species_name = get_taxa_name()

    retain = check_completeness()
    create_good_set(species_name, retain)

if __name__ == '__main__':
    main()
