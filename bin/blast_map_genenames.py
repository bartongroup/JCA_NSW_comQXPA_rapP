#!/usr/bin/env python

"""
Generates report on genes covered by a genomic blast hit
"""

import argparse
from pathlib import Path
import re
import sys

import json
import pandas as pd
from Bio import SearchIO

install_dir = Path(__file__).parent.parent
BAKTA_DIR = install_dir / 'data/full/annotations'

def get_match(string, match_re):

    """
    Searches a string with an re and returns the first match group

    Required params:
        string(str): query string 
        match_re(re): compiled re
    
    Returns:
        match(str): Contents of 1st match group
    """
    matches = re.search(match_re, string)

    if matches:
        gene_name = matches.group(1)
    else:
        gene_name = 'Unknown'

    return gene_name

def get_spanned_genes(assembly_accession, seq_accession, start, end):

    """
    Extracts names of genes spanned by hit

    Required params:
        assembly_accession(str): GCA accession 
        seq_accession(str): ID of sequence within assembly
        start(str): Start co-ordinate
        end(str): End co-ordinate
    """

    gff_file = BAKTA_DIR / assembly_accession / f"{assembly_accession}.gff3"
    gff_cols=('seqid', 'source', 'type', 'start', 'end', 'score', 'strand',
              'phase', 'attributes')
    gff = pd.read_csv(gff_file, sep="\t", comment='#', header=0, names=gff_cols, low_memory=False)

    gff = gff[(gff['seqid'].str.contains(seq_accession)) & (gff['type'] == 'CDS') &
              (gff['start'] >= int(start)) & (gff['end'] <= int(end))]

    gene_re = re.compile(r'gene=([a-zA-Z0-9]*)$')

    genes = gff.apply(lambda x: get_match(x['attributes'], gene_re), axis=1)

    if not genes.empty:
        genes = ','.join([x for x in genes])
        return genes

def parse_blast(record, mapping):

    """
    Parses a blast record
    
    Required parameters:
        record(pathlib.Path): Path to results file
        mapping(dict): Parsed seq id -> accession info mapping

    Returns:
        None
    """
    out_hits = []

    with record.open() as fh:
        with open('annotated.out','w', encoding='UTF-8') as out_fh:
            results = SearchIO.parse(fh, 'blast-tab')
            for result in results:

                # collect overall hit span 
                start_loc, end_loc = 0, 0

                for hit in result.hits:

                    seq_accession = hit.id.split('|')[2]
                    seq_info = mapping[seq_accession]
                    assembly_accession = seq_info['accession']
                    strain = seq_info['strain']

                    for hsp in hit.hsps:
                        # we don't care about no steenking strand...
                        start = min(hsp.hit_start, hsp.hit_end)
                        end = max(hsp.hit_start, hsp.hit_end)

                        start_loc = max(start, start_loc)
                        end_loc = max(end, end_loc)

                        pident = hsp.ident_pct
                        span = hsp.hit_span
                        mismatches = hsp.mismatch_num
                        gaps = hsp.gapopen_num
                        qstart = int(hsp.query_start) + 1
                        qend = hsp.query_end
                        hstart = int(hsp.hit_start) + 1
                        hend = hsp.hit_end
                        evalue = hsp.evalue
                        bitscore = hsp.bitscore

                    genes = get_spanned_genes(assembly_accession, seq_accession, start, end)

                    description = f"Assembly {assembly_accession}; Strain {strain};"
                    if genes:
                        description = f"{description}; Includes genes {genes}"
                    hit.description = description

                    out_line = f"{hit.query_id}\t{hit.id}\t{pident}\t{span}\t{mismatches}"
                    out_line = f"{out_line}\t{gaps}\t{qstart}\t{qend}\t{hstart}\t{hend}\t{evalue}\t{bitscore}\t{description}\n"
                    out_fh.write(out_line)


def main():
    """main method"""

    parser = argparse.ArgumentParser(
        prog = "blast_map_genenames.py",
        description="Reports genes which are spanned by a genomic blast hit"
    )
    parser.add_argument('-b', '--blast', required=True, help="Path to blast record")
    #parser.add_argument('-g', '--bakta', required = Treue, help='Path to bakta gff3 annotation directory')
    parser.add_argument('-m', '--mapping', required=True, help="Path to id mapping json file")

    args = parser.parse_args()
    result = Path(args.blast)

    with open(args.mapping, 'r', encoding='UTF-8') as fh:
        mapping = json.load(fh)


    parse_blast(result, mapping)

if __name__ == "__main__":
    main()