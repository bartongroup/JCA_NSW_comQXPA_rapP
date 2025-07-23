#!/usr/bin/env python

# There may be multiple copies of the operon present, so check for counts of annotated genes
#
# This is a bit tricky since comP or mutated genes are not necessarily going to have the correct
# names assigned

from pathlib import Path

from Bio import SeqIO
import pandas as pd

annot_dir = Path('data/full/annotations')
annots = annot_dir.glob("*/*.embl")

genes = {}
cds = {}

for genome in annots:
    accession = genome.stem

    gene_counts = {
        'comP': 0,
        'comQ': 0,
        'comA': 0,
        'comA2': 0,
        'comX': 0,
    }

    cds_counts = {
        'comP': 0,
        'comQ': 0,
        'comA': 0,
        'comA2': 0,
        'comX': 0,
    }

    for record in SeqIO.parse(genome, 'embl'):
        for feature in record.features:
            name = feature.qualifiers.get('gene')
            if name and name[0] in ['comP', 'comQ', 'comA', 'comA2', 'comX']:
                if feature.type == 'gene':
                    gene_counts[name[0]] += 1
                elif feature.type == 'CDS':
                    cds_counts[name[0]] += 1

    genes[accession] = gene_counts
    cds[accession] = cds_counts

gene_df = pd.DataFrame.from_dict(genes, orient="index")
cds_df = pd.DataFrame.from_dict(cds, orient='index')

df = gene_df.merge(cds_df, left_index = True, right_index = True, suffixes = ['_gene','_cds'])
df.to_csv('results/gene_counts.txt', sep="\t")