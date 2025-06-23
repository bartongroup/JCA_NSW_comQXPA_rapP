#!/usr/bin/env python

"""
Carries out series of pair-wise comparisons between Bacillus subtilis 168
sequences from subtiwiki and those isolated from each genome

Sequences are expected to be found in 'complete/fasta/com_proteins', generated
by `extract_sequences.py`. BS168 sequences are available in the
'reference_sequences' directory
"""

from subprocess import run
from argparse import ArgumentParser
from pathlib import Path
import re

from Bio import SeqIO
from Bio import SearchIO

import pandas as pd

parser = ArgumentParser(
    prog = __file__, 
        description="Identifies completion of sequences based upon comparison with 168 isolate"
    )
parser.add_argument('-m', '--metadata', required=True, help="Path to metadata file")
args = parser.parse_args()

reference_dir = Path('reference_sequences')
protein_dir   = Path('data/full/fasta/com_proteins')

metadata_df = pd.read_csv(args.metadata, sep="\t")

gene_results = dict()

for gene in ['comA', 'comP', 'comX', 'comQ']:
    results = list()

    query = reference_dir / Path(f'{gene}.fa')
    with open(query,'r', encoding='UTF-8') as fh:
        for record in SeqIO.parse(fh, format = 'fasta'):
            query_len = len(record.seq)
            break

    genomes = Path('data/full/fasta/genomes').glob('*fasta')

    for genome in genomes:
        # TBLASTN of known protein vs genome to see if it is present
        # The hits will be variable for i.e. comP which are not well conserved

        command = ['tblastn','-query', query, '-subject', genome, '-outfmt', '5', '-out', '/tmp/result']
        run(command, check=True)

        bl_result = SearchIO.read('/tmp/result', 'blast-xml')
        if len(bl_result.hits):
            genome = bl_result.hits[0].id.split('|')[1]

            strand = '+'
            if '-' in str(bl_result.hits[0].hsps[0].hit_frame):
                strand = '-'

            hit_res = {
                     'genome': genome,
                     'description': bl_result.hits[0].description.split(' ',1)[1],
                     'genome_query_len': query_len,
                     'genome_align_len': bl_result.hits[0].hsps[0].aln_span,
                     'genome_coverage': query_len/bl_result.hits[0].hsps[0].aln_span * 100,
                     'genome_identity': bl_result.hits[0].hsps[0].ident_num,
                     'genome_coords': f"{bl_result.hits[0].hsps[0].hit_start}-{bl_result.hits[0].hsps[0].hit_end}",
                     'genome_strand': strand
                 }
        else:
            hit_res = {'genome': genome}

        # BLASTP vs extracted sequence to make sure the annotation has found it correctly

        protein_file = protein_dir / gene / f"{genome}.fasta"
        if protein_file.exists():

            command = ['blastp', '-query', query, '-subject', protein_file, '-outfmt', '5', '-out', '/tmp/result']
            run(command, check = True)

            bl_result = SearchIO.read('/tmp/result', 'blast-xml')
            if len(bl_result.hits):

                hit_res['protein_align_len'] =  bl_result.hits[0].hsps[0].aln_span
                hit_res['protein_coverage'] =  query_len/bl_result.hits[0].hsps[0].aln_span * 100
                hit_res['protein_identity'] = bl_result.hits[0].hsps[0].ident_num

        # Flag cases where we seem to have better blast hits to genome than from identifiedo prediction
        hit_res['warning'] = False
        if ('protein_coverage' not in hit_res)  or \
             ('genome_coverage' not in hit_res) or \
             (hit_res['protein_coverage'] < hit_res['genome_coverage']) or \
             hit_res['genome_coverage'] < 90:
           hit_res['warning'] = True

        results.append(hit_res)

    results_df = pd.DataFrame(results)
    gene_results[gene] = results_df

merged_df = pd.merge(metadata_df, gene_results['comA'], suffixes=[None, '_comA'], left_on='Accession', right_on = 'genome', how = 'left')
merged_df.rename(
    columns = {'genome_query_len': 'genome_query_len_comA', 
                'genome_align_len': 'genome_align_len_comA', 
                'genome_coverage': 'genome_coverage_comA', 
                'genome_identity': 'genome_identity_comA', 
                'genome_coords': 'genome_coords_comA', 
                'protein_align_len': 'protein_align_len_comA', 
                'protein_coverage': 'protein_coverage_comA', 
                'protein_identity': 'protein_identity_comA', 
                'warning': 'warning_comA'
                }, inplace = True)

merged_df = pd.merge(merged_df,gene_results['comP'], suffixes=[None,'_comP'], on='genome', how = 'left')
merged_df = pd.merge(merged_df,gene_results['comX'], suffixes=[None, '_comX'], on='genome', how = 'left')
merged_df = pd.merge(merged_df,gene_results['comQ'], suffixes=[None, '_comQ'], on='genome', how = 'left')
merged_df.drop(['genome','description', 'description_comP','description_comX', 'description_comQ'], axis = 1, inplace=True)

with pd.ExcelWriter('results/protein_coverage.xlsx') as writer:
    merged_df.to_excel(writer, sheet_name = 'protein coverage', index = False)

