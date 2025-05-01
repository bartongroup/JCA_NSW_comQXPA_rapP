#!/usr/bin/env python

"""
Carries out series of pair-wise comparisons between Bacillus subtilis 168
sequences from subtiwiki and those isolated from each genome

Sequences are expected to be found in 'complete/fasta/com_proteins', generated
by `extract_sequences.py`. BS168 sequences are available in the
'reference_sequences' directory
"""

from pathlib import Path
from Bio import SeqIO
from Bio import SearchIO
from subprocess import run

import pandas as pd
import re

reference_dir = Path('reference_sequences')
protein_dir   = Path('complete/fasta/com_proteins')
work_dir = Path('protein_check')

work_dir.mkdir(exist_ok = True)

gene_results = dict()

for gene in ['comA', 'comP', 'comX', 'comQ']:
    results = list()

    query = reference_dir / Path(f'{gene}.fa')
    with open(query,'r') as fh:
        for record in SeqIO.parse(fh, format = 'fasta'):
            query_len = len(record.seq)
            break

    with open(protein_dir / Path(f'{gene}.faa')) as fh:
        for record in SeqIO.parse(fh, format = 'fasta'):
            subject = work_dir / 'subject.fa'
            result = work_dir / 'result.out'
            subject_len = len(record.seq)

            with open(subject, 'w') as out_fh:
                SeqIO.write([record], out_fh, format = 'fasta')

            command = ['blastp', '-query', query, '-subject', subject, '-outfmt', '5', '-out', result]
            run(command, check = True)

            bl_result = SearchIO.read(result, 'blast-xml')
            if len(bl_result.hits):

                genome = re.sub(r'_[A-Za-z0-9]+$', '', bl_result.hits[0].id)
                ident = bl_result.hits[0].hsps[0].ident_num

                hit_res = {
                    'genome': genome,
                    'description': bl_result.hits[0].description.split(' ',1)[1],
                    'query_len': query_len,
                    'subject_len': subject_len,
                    'coverage': (subject_len / query_len) * 100,
                    'ident': (ident / subject_len) * 100
                }
            else:
                hit_res = {'genome': genome}
            results.append(hit_res)

            subject.unlink()
            result.unlink()


    results_df = pd.DataFrame(results)
    gene_results[gene] = results_df

work_dir.unlink()

merged_df = pd.merge(gene_results['comA'],gene_results['comP'], suffixes=['_comA','_comP'], on='genome', how = 'outer')
merged_df = pd.merge(merged_df,gene_results['comX'], suffixes=[None, '_comX'], on='genome', how = 'outer')
merged_df = pd.merge(merged_df,gene_results['comQ'], suffixes=[None, '_comQ'], on='genome', how = 'outer')

with pd.ExcelWriter('protein_coverage.xlsx') as writer:
    merged_df.to_excel(writer, sheet_name = 'protein coverage', index = False)

