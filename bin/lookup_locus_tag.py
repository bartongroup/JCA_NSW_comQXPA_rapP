#!/bin/env python

from argparse import ArgumentParser
from pathlib import Path
import gzip
import json
import re
import sys

from Bio import SeqIO

def parse_data(genome_dir, serialised_data):

    """
    Parses genome records to extract necessary accessions

    Requred arguments:
        genome_dir(Path): Path to genome directory
    
    Returns:
        accessions(dict): Dict of data keyed on locus tag
    """
    if not genome_dir.exists():
        raise FileNotFoundError(f"Directory '{genome_dir}' not found")

    compressed = True
    genomes = list(genome_dir.glob('*gz'))

    if len(genomes) == 0:
        genomes = list(genome_dir.glob('*embl'))
        compressed = False

    if len(genomes) == 0:
        raise RuntimeError(f"No genomes (.embl/.gz) found in '{genome_dir}'")

    data= {}
    for genome in genomes:

        genome_data = {}
        genome_data['genes'] = {}
        locus_tag = ""
        short_tag = ""

        if compressed:
            fh = gzip.open(genome, 'rt')
        else:
            fh = open(genome, encoding = 'UTF-8')

        for record in SeqIO.parse(fh, format = 'embl'):
            genome_accession = genome.name
            sequence_id = record.id
            genome_accession = re.sub('.embl(.gz)?', '', genome_accession)

            for feature in record.features:
                if feature.type == 'source':

                    genome_data['organism'] = feature.qualifiers['organism'][0]
                    if 'strain' in feature.qualifiers:
                        genome_data['strain'] = feature.qualifiers['strain']

                elif feature.type == 'CDS':

                    locus_tag = feature.qualifiers['locus_tag'][0]
                    short_tag = locus_tag.split('_')[0]

                    if 'locus_tag' not in genome_data:
                        genome_data['locus_tag'] = locus_tag.split('_')[0]
                        genome_data['accession'] = genome_accession

                    genome_data['genes'][locus_tag] = {}
                    for qual in ('gene', 'product', 'protein_id'):
                        genome_data['genes'][locus_tag][qual] = \
                            feature.qualifiers[qual][0] if qual in feature.qualifiers else None
                    genome_data['genes'][locus_tag]['sequence_id'] = sequence_id
                    genome_data['genes'][locus_tag]['coordinates'] = f"{feature.location.start} - {feature.location.end}"
                    genome_data['genes'][locus_tag]['strand'] = feature.location.strand

        fh.close()

        if short_tag:
            data[short_tag] = genome_data

    with gzip.open(serialised_data, 'w') as out_fh:
        out_fh.write(json.dumps(data).encode('UTF-8'))

    return data

def print_details(ids, data):
    """
    Outputs summary of protein for each identifier in input

    Required params
        ids(string): (possibly multi-line) containing ids
        data(dict): Dictionary keyed on locus tax prefix produced by parse_data()

    Returns: 
        None
    """

    id_list = ids.split()

    for locus_id in id_list:

        tag = locus_id.split('_')[0]
        locus_data = data[tag]['genes'][locus_id]

        fields = (locus_id, data[tag]['accession'], data[tag]['organism'], 
                  f"gene: {locus_data['gene']}", locus_data['product'], locus_data['sequence_id'], 
                  locus_data['coordinates'], str(locus_data['strand']))

        print("\t".join(fields))

def main():

    """
    Main method...
    """

    parser = ArgumentParser(
        prog = __file__,
            description = """
            Looks up genome/gene associated with a locus tagged identifier i.e.
            HLBEOP_22245. If preparsed data is available (locus_tag_data.json)
            this will be used, otherwise the genomes in the directory specified
            by '--genomes' will be read, and the results stored for reuse""" )

    parser.add_argument('-i', '--input', help="Idenfiers(s) to resolve")
    parser.add_argument('-j', '--json', help="Path to genome_locus_tags.json.gz file", 
        default = '/cluster/mmb/common/Bsubtilis_complete_genomes/validated/genome_locus_tags.json.gz')
    args = parser.parse_args()

    serialised_data = Path(args.json)
    genome_dir = serialised_data.parent / 'genomes'

    if serialised_data.exists():

        with gzip.open(serialised_data, 'rb') as fh:
            json_data = fh.read()
        data = json.loads(json_data)

    else:
        data = parse_data(genome_dir, serialised_data)

    if args.input:
        input = args.input
    else:
        input = sys.stdin.read()

    print_details(input, data)

if __name__ == "__main__":
    main()
