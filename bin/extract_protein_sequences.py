#!/usr/bin/env python

'''
Extracts genes of interest from reannotated genomes

Script specific to comP - not intended to be generic
'''

import argparse
import re
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

NCBI_TAXID = '1423'

def get_qualifier(name, feature):

    """
    Extracts a qualifier value from a feature for a given name

    Required params:
        name(str): qualifier name
        feature(Bio::SeqFeature): feature of interest
    """

    if feature.qualifiers.get(name):
        return feature.qualifiers.get(name)[0]

    return None


def extract_feature_details(feature):

    """
    Pulls required fields from feature object

    Required parameters: 
        feature(Bio::SeqFeature): feature of interest
    
    returns:
        feature_info(dict): Dictionary containing required fields
    """

    feature_info = {
        'location':  f"{feature.location.start}-{feature.location.end}"
    }

    translation = get_qualifier('translation', feature)
    feature_info['translation'] = translation
    if translation is not None:
        feature_info['cds_length'] = len(translation)

    feature_info['gene_id'] = get_qualifier('gene', feature)
    feature_info['product'] = get_qualifier('product', feature)

    feature_info['pseudogene'] = get_qualifier('pseudogene', feature)

    if feature.qualifiers.get('note') and feature_info['pseudogene']:
        feature_info['note'] = feature.qualifiers.get('note')[-1]
    else:
        feature_info['note'] = None

    feature_info['locus_tag'] = get_qualifier('locus_tag', feature)

    return feature_info

def extract_cds_details(cds_info, record, accession, gene_name):

    """
    Extracts details from CDS features, including creating sequence
    object from translation

    Required fields:
        cds_info(dict): Dict of details of CDS
        record(SeqRecord): Parent record for CDS feature
        accession(str): Genome accession

    Returns:
        cds_info(dict): Updated dict
        protein(Bio::SeqRecord): Sequence of protein translation
    """

    protein = None

    cds_info['record_id'] = record.id
    cds_info['genomic_context'] = 'chromosome'
    if 'plasmid' in record.description:
        cds_info['genomic_context'] = 'plasmid'
    gene_id = cds_info.get("gene_id")
    if gene_id is None:
        gene_id = gene_name
    gene_id = f'{accession}_{gene_id}'

    if cds_info.get('translation'):

        protein = SeqRecord(
            Seq(cds_info['translation']),
                id = gene_id,
                description = cds_info.get('product'),
                annotations={'molecule_type': "protein"}
        )

    return(cds_info, protein)


def get_gene_sequences(genome):

    """
    Parses the downloaded genome record to extract the required sequences

    required parameters:
        metadata(pd.DataFrame): DataFrame mapping accessions to strain info
        genome(libpath.Path): Path to record

    returns:
        accession(str): accession of record
        gene_info(dict): metadata dictionary keyed on gene name
        cds_seqs(dict): dict of Bio::Seq objects keyed on gene name
    """

    gene_info = {}
    cds_seqs = {}
    prot_seqs = {}

    accession = str(genome.stem).replace('.embl','')
    with open(genome, 'rt', encoding='UTF-8') as fh:

        for record in SeqIO.parse(fh, format = 'embl'):

        #comP and comQ are not reliably annotated, however comX is present in
        #every annotation so use this to identify the location of the operon

            for index, feature in enumerate(record.features):

                if feature.type == 'gene' and 'gene' in feature.qualifiers.keys():
                    if 'comX' in feature.qualifiers.get('gene')[0]:
                        comX_gene_index = index
                        cds_indices = {}

                        # The operon is generally on the -ve strand, but there are
                        # some cases where it is positive

                        # Since we have consistent annotation formats, we should
                        # be able to rely on the surrounding features being the
                        # comX CDS, and the comP and comQ genes and CDSs, which
                        # helps given some sequences are fairly divergent

                        if feature.location.strand== -1:

                            gene_info['strand']='-'
                            cds_indices['comX'] = comX_gene_index + 1
                            cds_indices['comP'] = comX_gene_index - 1
                            cds_indices['comQ'] = comX_gene_index + 3
                            cds_indices['comA'] = comX_gene_index - 3
                        else:

                            gene_info['strand']='+'
                            cds_indices['comX'] = comX_gene_index + 1
                            cds_indices['comP'] = comX_gene_index + 3
                            cds_indices['comQ'] = comX_gene_index - 1
                            cds_indices['comA'] = comX_gene_index + 5

                        for index, feature in cds_indices.items():

                            cds_info = extract_feature_details(record.features[feature])

                            # comP frameshifts can not only cause a truncation, but also a second product from the 3' region,
                            # which screws up comA positioning, as do transposase insertions within comP

                            # These checks work with all available genomes at the time of writing but changes may be required
                            # in the light of novel sequences being produced
                            if index == 'comA':
                                if not cds_info['product']:
                                    print(f"{record.id}: No product available" )
                                elif 'histidine kinase' in cds_info['product'] or 'hypothetical protein' in cds_info['product']:
                                    if gene_info['strand']=='-':
                                        cds_info = extract_feature_details(record.features[comX_gene_index-5])
                                    else:
                                        cds_info = extract_feature_details(record.features[comX_gene_index+7])
                                elif 'transposase' in cds_info['product']:
                                    if gene_info['strand']=='-':
                                        cds_info = extract_feature_details(record.features[comX_gene_index-7])
                                    else:
                                        cds_info = extract_feature_details(record.features[comX_gene_index+9])

                            gene_info[index], prot_seqs[index] = extract_cds_details(cds_info, record, accession, index)

                # rapP does not seem to be annotated with a gene name by bakta, but does have
                # a product attached to CDS features
                elif feature.type == 'CDS' and 'product' in feature.qualifiers.keys() and \
                    feature.qualifiers.get('product')[0] == 'Rap phosphatase':

                    cds_info = extract_feature_details(feature)
                    gene_info['rapP'], prot_seqs['rapP'] = extract_cds_details(cds_info, record, accession, 'rapP')

    return(accession, gene_info, cds_seqs, prot_seqs)

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
        df = df[['accession','record_id','genomic_context', 'strain','gene_id','locus_tag','cds_length','product','strand','location','pseudogene','note']]
        dfs[gene] = df

    with pd.ExcelWriter(f'{outpath}/protein_status.xlsx') as writer:
        for gene, gene_df in dfs.items():
            gene_df.to_excel(writer, sheet_name = gene, index = False)

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

        accession, gene_info, cds_seqs, prot_seqs = get_gene_sequences(genome)

        all_cds_seqs.append(cds_seqs)
        all_prot_seqs.append(prot_seqs)
        all_gene_info[accession] = gene_info

    write_seqs(all_prot_seqs, Path(f'{args.datadir}/fasta/target_proteins'))
    summarise_genes(metadata, Path(args.datadir), all_gene_info)

if __name__ == "__main__":
    main()
