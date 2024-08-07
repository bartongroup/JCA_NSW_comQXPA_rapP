#!/bin/env python

"""
Extracts operon sequence as far as is possible based on annotated feature
co-ordinates in the 'complete/gene_summary.xlsx' sheet. We are looking
for the range from the start of comQ to the end of comA

Adds additional column to spreadsheet of genomes including the length of the
isolated operon
"""

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from tqdm import tqdm

GENOME_DIR = Path('refined/genomes')
OPERON_DIR = Path('operon')

def main():

    comQ_info = pd.read_excel('complete/gene_summary.xlsx', sheet_name = 'comQ')
    comA_info = pd.read_excel('complete/gene_summary.xlsx', sheet_name = 'comA')

    # beware fuzzy co-ordinates...
    comQ_info.replace(r'[<>]','',regex = True, inplace= True)
    comA_info.replace(r'[<>]','',regex = True, inplace= True)

    lengths = list()
    output_records = list()

    for genome in tqdm(GENOME_DIR.glob('*.embl')):

        record_id = comQ_info[comQ_info['accession'] == genome.stem]['record_id'].values[0]
        strand    = comQ_info[comQ_info['accession'] == genome.stem]['strand'].values[0]
        comQ_loc  = comQ_info[comQ_info['accession'] == genome.stem]['location'].values[0].split('-')
        comA_loc  = comA_info[comA_info['accession'] == genome.stem]['location'].values[0].split('-')

        if strand == '+':
            start = int(comQ_loc[0])
            end   = int(comA_loc[1])
        else:
            start = int(comA_loc[0])
            end   = int(comQ_loc[1])

        span = end - start

        with open(genome) as fh:
            for record in SeqIO.parse(fh, 'embl'):
                if record.id == record_id:
                    operon_seq = record.seq[start:end]
                    if strand == '-':
                        operon_seq = operon_seq.reverse_complement()
                    output_record = SeqRecord(id = genome.stem, description='comQ_comX_comP_comA operon', seq = operon_seq)
                    output_records.append(output_record)
                else:
                    next

        lengths.append({'accession': genome.stem, 'operon length': span, 'reference difference': span - 4188})

    with open(OPERON_DIR / "operons.fasta", 'w') as out_fh:
        SeqIO.write(output_records, out_fh, 'fasta')
    
    genome_info = pd.read_excel('Bacillus_subtilis_complete_genomes_25-06-2024.xlsx', sheet_name = 'Complete genomes')
    length_df = pd.DataFrame(lengths)
    genome_info = pd.merge(genome_info, length_df, how = 'right', left_on = 'accession', right_on = 'accession')

    with pd.ExcelWriter('refined/refined_genomes.xlsx') as writer:
        genome_info.to_excel(writer, sheet_name = 'Refined genomes', index = False)

if __name__ == "__main__":
    main()