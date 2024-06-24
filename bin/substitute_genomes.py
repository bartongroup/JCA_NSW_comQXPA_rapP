#!/bin/env python

"""
Improved versions of some of the genomes have been shared by the group
responsible for them, so where possible these need to be replaced in the 
database
"""

import pandas as pd
from Bio import SeqIO
from pathlib import Path
from shutil import copyfile
from tqdm import tqdm

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
    danish_dir =  Path(__file__).parents[1] / Path('Danish_sequences')

    excel_files = list(Path('.').glob('[!.]*.xlsx'))
    genome_df = pd.read_excel(excel_files[-1], dtype={'isolate': str})

    genome_df['munged_isolate']=genome_df['isolate'].apply(lambda x: x.replace(' ',''))
    match_count = 0
    novel_genomes = list()

    for genome in Path('Danish_sequences').glob('Bacillus_subtilis_*.fasta'):
        isolate = genome.name.replace('Bacillus_subtilis_','').replace('.fasta','')

        matched = genome_df[genome_df.munged_isolate.str.contains(f"{isolate}$", regex = True)]

        if not matched.empty:
            match_count+=1
            accession = matched['accession'].values[0]
            ena_isolate = matched['isolate'].values[0]
            name = genome.name
            stem = genome.stem

            replace_genome(fasta_dir, blast_dir, accession, name)
            genome_df.loc[genome_df["isolate"]==ena_isolate, 'accession'] = stem
        else: 
            novel_genomes.append(name)

    novel_genomes = list(set(novel_genomes))

    for genome in novel_genomes:

        stem = Path(genome).stem
        scaffolds = 0
        with open(danish_dir / Path(genome), 'r') as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                scaffolds = scaffolds + 1

        copyfile( Path('Danish_sequences', genome), fasta_dir / Path(genome).name )

        # add details to genomes dataframe
        genome_dat = {
            'accession': stem,
            'isolate':   stem,
            'source':    None,
            'scaffolds': scaffolds
        }
        new_df = pd.DataFrame(genome_dat, index = [0])
        genome_df = pd.concat([genome_df, new_df])

    # Blast index the full dataset again...
    blast_index(fasta_dir, blast_dir, 'Bacillus subtilis')
    genome_df.drop(['munged_isolate'], inplace = True, axis = 1)
    genome_df.to_excel(excel_files[-1], sheet_name = 'Complete genomes', index = False, header = True)
    genome_df.drop(['source', 'scaffolds'], inplace = True, axis = 1)
    genome_df.to_csv('strains.txt', sep="\t", header = True, index = False)

if __name__ == "__main__":
    main()