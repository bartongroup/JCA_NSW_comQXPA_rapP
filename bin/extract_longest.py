#!/bin/env python

"""
Attempts to separate the chromosomal sequence from any plasmids or unintegrated
ribosomal clusters  to provide a more meaningful ANI

Works from refined genomes only and outputs fasta records of longest sequences

NB - turns out mutliple fragments are not leading to low ANI scores.
"""

from Bio import SeqIO
from pathlib import Path

GENOME_DIR = Path('refined/genomes')
LONG_RECORD_DIR = Path('refined/long_sequences')

def main():

    LONG_RECORD_DIR.mkdir(exist_ok = True)

    genomes = GENOME_DIR.glob('*.embl')

    for genome in genomes:
        records = SeqIO.to_dict(SeqIO.parse(genome, 'embl'))

        if len(records) > 1:
            longest = None
            for record in records.values():

                if (longest is None) or (len(record.seq) > len(longest.seq)):
                    longest = record

            out_record = longest
            if len(longest.seq) < 4000000:
                print(f"WARNING: {genome.stem}: {len(out_record.seq)}")
        else:
            out_record = list(records.values())[0]
        
        with open(LONG_RECORD_DIR / f"{genome.stem}.fasta", 'w' ) as out_fh:
            SeqIO.write([out_record], out_fh, 'fasta')


if __name__ == "__main__":
    main()