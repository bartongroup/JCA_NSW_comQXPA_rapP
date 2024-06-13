#!/bin/env python

"""
Improved versions of some of the genomes have been shared by the group
responsible for them, so where possible these need to be replaced in the 
database
"""

import pandas as pd
from pathlib import Path
from shutil import copyfile

from build_genome_dbs import blast_index

def replace_genome(fasta_dir, blast_dir, accession, name):

    """
    Replaces an ENA genome in the 'fasta' directory with the improved one
    and also replaces the blast database

    required parameters:
        accession(str): Accession of genome to remove
        name(str): Genome to replace with

    returns:
        None
    """

    blast_suffixes = ['.ndb', '.nhr', '.nin', '.njs', '.not', '.nsq', '.ntf', '.nto']

    new_genome = Path('Danish_sequences', name)
    fasta_file = fasta_dir / ( accession + '.fasta')

    if fasta_file.exists():
        fasta_file.unlink()

    for suffix in blast_suffixes:

        blast_file = blast_dir / (accession + suffix)

        if blast_file.exists():
            blast_file.unlink()
        else:
            print(f'{blast_file} not found')

    copyfile(new_genome, fasta_dir / name)


def main():

    fasta_dir  = Path(__file__).parents[1] / Path('fasta')
    blast_dir  = Path(__file__).parents[1] / Path('blast_db')

    genome_df = pd.read_excel('Bacillus_subtilis_complete_genomes_03-06-2024.xlsx')
    genome_df['isolate']=genome_df['isolate'].apply(lambda x: x.replace(' ',''))

    match_count = 0
    novel_genomes = list()

    for genome in Path('Danish_sequences').glob('Bacillus_subtilis_*.fasta'):
        isolate = genome.name.replace('Bacillus_subtilis_','').replace('.fasta','')

        matched = genome_df[genome_df.isolate.str.contains(f"{isolate}$", regex = True)]

        if not matched.empty:
            match_count+=1
            accession = matched['accession'].values[0]
            ena_isolate = matched['isolate'].values[0]
            name = genome.name

            replace_genome(fasta_dir, blast_dir, accession, name)

        else: 

            novel_genomes.append(name)
    
    for genome in novel_genomes:

        copyfile( Path('Danish_sequences', genome), fasta_dir / Path(genome).name )
    
    # Blast index the full dataset again...
    blast_index(fasta_dir, blast_dir, 'Bacillus subtilis')

if __name__ == "__main__":
    main()