#!/usr/bin/env python

'''
Extracts genes of interest from reannotated genomes

Script specific to comP - not intended to be generic
'''

import argparse
from pathlib import Path

from Bio import SeqIO

from common import get_target_sequences

def main():

    """
    Main method
    """

    parser = argparse.ArgumentParser(
        prog = __file__,
        description="Extracts target sequences from outgroup genome"
    )
    parser.add_argument('-d', '--datadir', required=True, help="Path to data directory")
    args = parser.parse_args()

    genomes = list(Path(f"{args.datadir}").glob('*.embl'))

    accession, gene_info, cds_seqs, prot_seqs = get_target_sequences(genomes[0])

    target_path = Path(f'{args.datadir}/target_proteins')
    target_path.mkdir(exist_ok = True)
    for protein_name, seq_record in prot_seqs.items():

        with open(Path(target_path / f"{protein_name}.fasta"), 'w', encoding='UTF-8') as fh:
            SeqIO.write([seq_record], fh, 'fasta')

if __name__ == "__main__":
    main()
