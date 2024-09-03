#!/bin/env python

"""
A final round of genome refinement to remove sequences identified by 
GTDBTK as being something other than B.subtilis

Required arguments: None
"""

import pandas as pd
from pathlib import Path
from os import symlink


GENOME_DIR = Path('refined/genomes')
FASTA_DIR = Path('refined/fasta/genomes')
CDS_DIR = Path('refined/fasta/cds')
PROTEIN_DIR = Path('refined/fasta/protein')

FINAL_GENOME_DIR = Path('final/genomes')
FINAL_FASTA_DIR = Path('final/fasta/genomes')
FINAL_CDS_DIR = Path('final/fasta/cds')
FINAL_PROTEIN_DIR = Path('final/fasta/protein')

def copy_genomes(accession):

    """
    Creates  symlink of necessary files in 'final' directory 
    from refined data set

    Required parameters:
        accession(str): Genome identifier
    
    Returns:
        None
    """

    genome = "../.." / GENOME_DIR / f"{accession}.embl"
    fasta = "../../.." / FASTA_DIR / f"{accession}.fasta"
    cds = "../../.." / CDS_DIR / Path(f"{accession}.ffn")
    protein = "../../.." /PROTEIN_DIR / Path(f"{accession}.ffn")

    symlink(genome, FINAL_GENOME_DIR / f"{accession}.embl")
    symlink(cds, FINAL_CDS_DIR / f"{accession}.ffn")
    symlink(protein, FINAL_PROTEIN_DIR / f"{accession}.ffn")
    symlink(fasta, FINAL_FASTA_DIR / f"{accession}.fasta")

def main():

    FINAL_GENOME_DIR.mkdir(exist_ok = True, parents = True)
    FINAL_FASTA_DIR.mkdir(exist_ok = True, parents = True)
    FINAL_CDS_DIR.mkdir(exist_ok = True, parents = True)
    FINAL_PROTEIN_DIR.mkdir(exist_ok = True, parents = True)

    outlier_df = pd.read_csv('results/gtdbtk/outlier_classifications.txt', sep="\t")
    outliers = outlier_df['accession'].values

    genomes = list(GENOME_DIR.glob('*.embl'))
    genomes = [x.name.replace('.embl','') for x in genomes]
    genomes = [x for x in genomes if x not in outliers]

    for accession in genomes:
        copy_genomes(accession)



if __name__ == "__main__":
    main()