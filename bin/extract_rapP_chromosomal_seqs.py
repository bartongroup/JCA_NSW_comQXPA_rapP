#!/usr/bin/env python3

"""
Extracts rapP carrying chromosome seqquences from annotated genomes
"""

from argparse import ArgumentParser
from pathlib import Path
from subprocess import run, CalledProcessError
from Bio import SeqIO
import sys
import pandas as pd

def check_accession(accession, chromosomal_metadata):
    """
    Check a given accession to identify the seqence of any rapP carrying chromosome sequences

    Required parameters:
        accession: GCA accession string
        chromosomal_metadata: DataFrame containing metadata on chromosomal rapP sequences

    Returns:
        info(dict): information on plasmid sequence and identity to pBS32 
    """
    print("Checking accession:", accession)

    genome_embl = Path(f'data/refined/annotations/{accession}/{accession}.embl')

    if genome_embl.exists():
        records = SeqIO.parse(genome_embl, format='embl')
        for record in records:
            if not 'plasmid' in record.description.lower():
                features = [f for f in record.features if f.type == 'CDS']
                for feature in features:
                    if 'rap phosphatase' in feature.qualifiers.get('product', [''])[0].lower():
                        print('Found rapP gene...extracting sequence') 

                        output_fasta = Path(f'data/refined/rapP_chromosomes/{accession}_rapP_chromosome.fasta')
                        output_embl = Path(f'data/refined/rapP_chromosomes/{accession}_rapP_chromosome.embl')

                        with open(output_fasta, 'w', encoding='UTF-8') as out_fh:
                            SeqIO.write(record, out_fh, format='fasta')

                        with open(output_embl, 'w', encoding='UTF-8') as out_fh:
                            SeqIO.write(record, out_fh, format='embl')  

                        print('Extracting rapP locus region...')
                        start, end = chromosomal_metadata.loc[chromosomal_metadata['Accession'] == accession, 'location_rapP'].values[0].split('-')
                        start = int(start) - 10000
                        end = int(end) + 10000

                        record = record[start:end]

                        locus_fasta = Path(f'data/refined/rapP_loci/{accession}_rapP_region.fasta')
                        locus_embl = Path(f'data/refined/rapP_loci/{accession}_rapP_region.embl')

                        with open(locus_fasta, 'w', encoding='UTF-8') as out_fh:
                            SeqIO.write(record, out_fh, format='fasta')

                        with open(locus_embl, 'w', encoding='UTF-8') as out_fh:
                            SeqIO.write(record, out_fh, format='embl')  

def main():
    """
    Main method
    """

    parser = ArgumentParser(
        prog = __file__,
        description="Extract rapP carrying chromosome sequences from annotated genomes"
    )

    parser.add_argument('-m', '--metadata', required=True, help="Path to metadata excel file")
    args = parser.parse_args()

    metadata = pd.read_excel(args.metadata, header=1)
    rapP_chromosomes = metadata[metadata['genomic_context_rapP'] == 'chromosome']
    chromosomal_accessions = rapP_chromosomes['Accession'].tolist()
    chromosomal_metadata = rapP_chromosomes[['Accession', 'strain_rapP', 'locus_tag_rapP', 'strand_rapP', 'location_rapP', 'cds_length_rapP']]

    chromosomal_metadata.to_csv('rapP_context/rapP_chromosomal_metadata.txt', sep="\t", index=False)

    results = []
    for accession in chromosomal_accessions:
        info = check_accession(accession, chromosomal_metadata)
        if info:
            results.append(info)

if __name__ == '__main__':
    main()