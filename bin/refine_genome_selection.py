#!/bin/env python

"""
The selection of genomes has been refined by NSW to remove duplicate
representation of isolates and those involved in genome reduction experiments.

This script subsets the complete genome set based upon these refinements
"""

import pandas as pd
from pathlib import Path
from shutil import copy

GENOME_DIR = Path('complete/genomes')
FASTA_DIR = Path('complete/fasta/genomes')
CDS_DIR = Path('complete/fasta/cds')
PROTEIN_DIR = Path('complete/fasta/protein')

REFINED_GENOME_DIR = Path('refined/genomes')
REFINED_FASTA_DIR = Path('refined/fasta/genomes')
REFINED_CDS_DIR = Path('refined/fasta/cds')
REFINED_PROTEIN_DIR = Path('refined/fasta/protein')

def refine_genomes(accession):

    """
    Creates copy of necessary files in 'refined' directory

    Required parameters:
        accession(str): Genome identifier
    
    Returns:
        None
    """

    genome = GENOME_DIR / Path(f'{accession}.embl')
    fasta = FASTA_DIR / Path(f'{accession}.fasta')
    cds = CDS_DIR / Path(f'{accession}.ffn')
    protein = PROTEIN_DIR / Path(f'{accession}.ffn')

    copy(genome, REFINED_GENOME_DIR)
    copy(fasta, REFINED_FASTA_DIR)
    copy(cds, REFINED_CDS_DIR)
    copy(protein, REFINED_PROTEIN_DIR)


def main():

    REFINED_GENOME_DIR.mkdir(exist_ok = True, parents = True)
    REFINED_FASTA_DIR.mkdir(exist_ok = True, parents = True)
    REFINED_CDS_DIR.mkdir(exist_ok = True, parents = True)
    REFINED_PROTEIN_DIR.mkdir(exist_ok = True, parents = True)

    genome_df = pd.read_excel('ComP_gene_summary_July.xlsx', sheet_name = 'finalComP', )
    genome_df['accession'].apply(refine_genomes)

    protein_coverage_df = pd.read_excel('protein_coverage.xlsx', 'protein coverage')
    protein_coverage_df = protein_coverage_df[protein_coverage_df['genome'].isin(genome_df['accession'])]

    with pd.ExcelWriter('protein_coverage_refined.xlsx') as writer:
        protein_coverage_df.to_excel(writer, sheet_name = 'protein coverage', index = False)

if __name__ == "__main__":
    main()