#!/usr/bin/env python

"""
Improved versions of some of the genomes have been shared by the group
responsible for them, so where possible these need to be replaced in the 
database
"""

# TODO: This is a temporary script required to add the Danish sequences to the dataset
# but needs to be removed once these are in the public databases and can be picked up 
# by the initial automated database build

from pathlib import Path
from shutil import copyfile

import pandas as pd
from Bio import SeqIO

from common import blast_index

NCBI_TAXID = 1423

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

        blast_file = blast_dir / 'genomes' / (accession + suffix)

        if blast_file.exists():
            blast_file.unlink()
        else:
            print(f'{blast_file} not found')

    copyfile(new_genome, fasta_file)


def main():

    """ main function """

    fasta_dir  = Path(__file__).parents[1] / Path('data/full/fasta/genomes')
    blast_dir  = Path(__file__).parents[1] / Path('data/full/blast_db/')
    danish_dir =  Path(__file__).parents[1] / Path('Danish_sequences')

    metadata_files = list(Path('data/full/').glob('Bacillus_subtilis_complete_genome*.txt'))
    metadata_file = metadata_files[-1]
    genome_df = pd.read_csv(metadata_file, dtype={'Strain': 'str'}, sep="\t")
    # don't know why this isn't being imported as a str given the dtype in the above line
    genome_df['Strain'] = genome_df['Strain'].apply(str)

    genome_df['munged_isolate']=genome_df['Strain'].apply(lambda x: x.replace(' ',''))
    match_count = 0
    novel_genomes = []

    for genome in danish_dir.glob('Bacillus_subtilis_*.fasta'):
        isolate = genome.name.replace('Bacillus_subtilis_','').replace('.fasta','')

        matched = genome_df[genome_df.munged_isolate.str.contains(f"{isolate}$", regex = True)]

        if not matched.empty:
            match_count+=1
            accession = matched['Accession'].values[0]
            ena_isolate = matched['Strain'].values[0]
            name = genome.name
            stem = genome.stem

            replace_genome(fasta_dir, blast_dir, accession, name)
            genome_df.loc[genome_df["Strain"]==ena_isolate, 'Danish ID'] = stem
        else:
            novel_genomes.append(name)

    novel_genomes = list(set(novel_genomes))

    for genome in novel_genomes:

        stem = Path(genome).stem
        scaffolds = 0
        with open(danish_dir / Path(genome), 'r', encoding='UTF-8') as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                scaffolds = scaffolds + 1

        copyfile( Path('Danish_sequences', genome), fasta_dir / Path(genome).name )

        # add details to genomes dataframe
        genome_dat = {
            'Accession': stem,
            'Danish ID': stem,
            'Strain':   stem,
            'Isolation source': None,
            'Scaffolds': scaffolds
        }
        new_df = pd.DataFrame(genome_dat, index = [0])
        genome_df = pd.concat([genome_df, new_df])

    # Blast index the full dataset again...
    blast_index(fasta_dir, blast_dir, species = 'Bacillus subtilis', ncbi_taxid = NCBI_TAXID, index_type = 'nucl')
    metadata_file = str(metadata_file).replace('.txt','_Danish_additions.txt')
    genome_df.drop(['munged_isolate'], inplace = True, axis = 1)
    genome_df.to_csv(metadata_file, sep="\t", index = False, header = True)

if __name__ == "__main__":
    main()