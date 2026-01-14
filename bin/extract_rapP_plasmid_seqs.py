#!/usr/bin/env python3

"""
Extracts rapP carrying plasmid sequences from annotated genomes

rapP_context/plasmid_isolates.txt contains a list of GCA accessions where rapP
has been identified on a plasmid 
"""

from argparse import ArgumentParser
from pathlib import Path
from subprocess import run, CalledProcessError
from Bio import SeqIO
import sys
import pandas as pd

def get_identity(seq1, seq2):
    """
    Carry out global alignment between two sequences using EMBOSS stretcher
    (Myers-Miller algorithm).
    The query sequence is first rotated to match the orientation of pBS32 
    
    Required parameters:
        seq1: Path to first sequence file
        seq2: Path to second sequence file
    
    Returns:
        ident(float): percentage identity between the two sequences
    """
    ref_start="ATGAGTAAATATTTCACAGCT"

    cmd = ['rotate', '-s', ref_start, '-o', f"/tmp/{seq1.stem}_rotated.fa", str(seq1)]

    try:
        rotate_result = run(cmd, capture_output=True, text=True, check=True)
    except CalledProcessError as e:
        print("Error running rotate command:", rotate_result.stderr)
        sys.exit(1)

    cmd = ['stretcher', 
           '-asequence', f"/tmp/{seq1.stem}_rotated.fa",
           '-bsequence', str(seq2), 
           '-gapopen', '10', 
           '-gapextend', '1', 
           '-outfile', f'/tmp/{seq1.stem}_vs_{seq2.stem}_needle.txt']

    try:
        stretcher_result = run(cmd, capture_output=True, text=True, check=True)
    except CalledProcessError as e:
        print("Error running stretcher command")
        print(stretcher_result.stdout)
        print(stretcher_result.stderr)
        sys.exit(1)

    with open(f'/tmp/{seq1.stem}_vs_{seq2.stem}_needle.txt') as fh:
        for line in fh:
            if line.startswith('# Identity:'):
                identity_line = line.strip()
                identity_line = identity_line.replace('# Identity: ', '').strip()
                print(identity_line)

                return identity_line

def check_accession(accession, plasmid_ref):
    """
    Check a given accession to identify the seqence of any rapP carrying plasmid
    and compare to pBS32

    Required parameters:
        accession: GCA accession string
        plasmid_ref: Path to plasmid reference sequence file

    Returns:
        info(dict): information on plasmid sequence and identity to pBS32 
    """
    print("Checking accession:", accession)

    genome_embl = Path(f'data/refined/annotations/{accession}/{accession}.embl')
    if genome_embl.exists():

        records = SeqIO.parse(genome_embl, format='embl')
        for record in records:
            if 'plasmid' in record.description.lower() or len(record.seq) < 100000:
                plasmid_length = len(record.seq)
                print('Have a plasmid sequence (length: ', len(record.seq), ')')
                features = [f for f in record.features if f.type == 'CDS']
                for feature in features:
                    if 'rap phosphatase' in feature.qualifiers.get('product', [''])[0].lower():
                        print('Found rapP gene on plasmid:')

                        output_fasta = Path(f'data/refined/rapP_plasmids/{accession}_rapP_plasmid.fasta')
                        with open(output_fasta, 'w', encoding='UTF-8') as out_fh:
                            SeqIO.write(record, out_fh, format='fasta')

                        output_embl = Path(f'data/refined/rapP_plasmids/{accession}_rapP_plasmid.embl')
                        with open(output_embl, 'w', encoding='UTF-8') as out_fh:
                            SeqIO.write(record, out_fh, format='embl')  

                        identity = get_identity(output_fasta, plasmid_ref)

                        return {
                            'accession': accession,
                            'plasmid_length': plasmid_length,
                            'pBS32_identity': identity
                        }
    return {
        'accession': None,
        'plasmid_length': None,
        'pBS32_identity': None
    } 
    

def main():
    """
    Main method
    """

    parser = ArgumentParser(
        prog = __file__,
        description="Extract rapP carrying plasmid sequences from annotated genomes, and determine identity to pBS32"
    )
    parser.add_argument('-m', '--metadata', required=True, help="Path to metadata excel file")
    parser.add_argument('-p', '--plasmid', required=True, help="Path to plasmid reference sequence file")
    args = parser.parse_args()

    plasmid_ref = Path(args.plasmid)

    metadata = pd.read_excel(args.metadata, header=1)
    rapP_plasmids = metadata[metadata['genomic_context_rapP'] == 'plasmid']

    plasmid_accessions = rapP_plasmids['Accession'].tolist()
    plasmid_metadata = rapP_plasmids[['Accession', 'strain_rapP', 'locus_tag_rapP', 'strand_rapP', 'location_rapP', 'cds_length_rapP']]

    results = []
    for accession in plasmid_accessions:
        info = check_accession(accession, plasmid_ref)
        if info:
            results.append(info)

    results_df = pd.DataFrame(results)

    plasmid_metadata = pd.merge(plasmid_metadata, results_df, left_on='Accession', right_on='accession', how='left')
    plasmid_metadata.to_csv('rapP_context/rapP_plasmid_identities.tsv', sep='\t', index=False)

if __name__ == '__main__':
    main()