#!/usr/bin/env python

'''
Extracts genes of interest from reannotated genomes

Script specific to comP - not intended to be generic
'''

import argparse
import re
from pathlib import Path

from Bio import SeqIO
import pandas as pd

from common import get_target_sequences

def write_seqs(seqs, db_path):

    """
    Writes fasta files of collected sequences

    Provided data is expected to be a dictionary keyed on accession, with values
    being a further dictionary keyed on gene name, with values being Bio.Seq
    objects

    Filenames and suffixes will be determined automatically 

    Required parameters:
        seqs(dict): Dictionary as described above
        db_path(pathlib.Path): Path to directory to write sequences
    """

    db_path.mkdir(exist_ok = True)

    # reorganise list of dicts into a dict of lists...
    seq_lists = {}

    for dic in seqs:
        for key in dic:
            if key in seq_lists:
                seq_lists[key].append(dic[key])
            else:
                seq_lists[key] = [dic[key]]

    for gene, seqs in seq_lists.items():
        Path(db_path / gene).mkdir(exist_ok=True)
        for seq_record in seqs:
            if seq_record is not None:
                genome = re.sub("_[A-Za-z]*$", "", seq_record.id)

                outfile = db_path / gene / f"{genome}.fasta"
                with open(outfile, 'w', encoding='UTF-8') as fh:
                    SeqIO.write([seq_record], fh, 'fasta')

def summarise_genes(strain_info, outpath, gene_info):

    """
    Produces summary spreadsheet on identified genes

    Required parameters:
        strain_info(pd.DataFrame): Mapping of accession to strain names
        outpath(str): Path to output directory
        gene_info(dict): dictionary keyed on accession
    
    returns:
        None
    """

    # reorganise list of dicts into a dict of lists...
    seq_lists = {}

    for acc in gene_info.keys():

        for key in gene_info[acc].keys():

            if key == 'strand':
                strand = gene_info[acc][key]
            else:
                dat = gene_info[acc][key]
                dat['accession'] = acc
                dat['strand'] = strand
                dat['strain'] = strain_info.loc[strain_info['Accession']==acc,'Title'].values[0]

                if key in seq_lists:
                    seq_lists[key].append(gene_info[acc][key])
                else:
                    seq_lists[key] = [gene_info[acc][key]]

    dfs = {}
    for gene,seq_list in seq_lists.items():
        df = pd.DataFrame(seq_list)
        df = df.rename(columns = lambda x: f"{x}_{gene}")
        df = df.rename({f"accession_{gene}": "accession"}, axis=1)
        df = df[['accession',f'record_id_{gene}',f'genomic_context_{gene}', f'strain_{gene}',
            f'gene_id_{gene}',f'locus_tag_{gene}',f'cds_length_{gene}',f'product_{gene}',f'strand_{gene}',
            f'location_{gene}',f'pseudogene_{gene}',f'note_{gene}']]
        dfs[gene] = df

    merged_df = pd.merge(dfs['comA'], dfs['comP'], on='accession', how = 'left')
    merged_df = pd.merge(merged_df, dfs['comX'], on='accession', how = 'left')
    merged_df = pd.merge(merged_df, dfs['comQ'], on='accession', how = 'left')
    merged_df = pd.merge(merged_df, dfs['rapP'], on='accession', how = 'left')

    merged_df.to_csv('data/full/protein_status.txt', sep="\t", index=False)

def main():

    """
    Main method
    """

    parser = argparse.ArgumentParser(
        prog = __file__,
        description="Downloads, fasta converts and optionally blast indexes genome records for organism"
    )
    parser.add_argument('-d', '--datadir', required=True, help="Path to data directory")
    parser.add_argument('-m', '--metadata', required=True, help="Path to metadata file")
    args = parser.parse_args()

    genomes = list(Path(f"{args.datadir}/annotations").glob('*/*.embl'))

    all_cds_seqs = list()
    all_prot_seqs  = list()
    all_gene_info = {}

    metadata = pd.read_csv(args.metadata, sep="\t")

    for genome in genomes:

        accession, gene_info, cds_seqs, prot_seqs = get_target_sequences(genome)

        all_cds_seqs.append(cds_seqs)
        all_prot_seqs.append(prot_seqs)
        all_gene_info[accession] = gene_info

    write_seqs(all_prot_seqs, Path(f'{args.datadir}/fasta/target_proteins'))
    summarise_genes(metadata, Path(args.datadir), all_gene_info)

if __name__ == "__main__":
    main()
