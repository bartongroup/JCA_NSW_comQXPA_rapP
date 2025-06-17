#!/bin/env python

from argparse import ArgumentParser
from pathlib import Path
import gzip
import json
import re

from Bio import SeqIO


def parse_data(genome_dir):

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
            genome_accession = re.sub('.embl(.gz)?', '', genome_accession)

            for feature in record.features:
                if feature.type == 'source':

                    genome_data['organism'] = feature.qualifiers['organism']
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

        fh.close()

        if short_tag:
            data[short_tag] = genome_data

    with gzip.open('genome_locus_tags.json.gz', 'w') as out_fh:
        out_fh.write(json.dumps(data).encode('UTF-8'))

    return data


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
    parser.add_argument('-g', '--genomes', help="Path to genomes directory")
    args = parser.parse_args()

    serialised_data = Path('genome_locus_tags.json.gz')

    if serialised_data.exists():

        with gzip.open('genome_locus_tags.json.gz', 'rb') as fh:
            json_data = fh.read()
        data = json.loads(json_data)

    else:
        if args.genomes:
            data = parse_data(Path(args.genomes))

        else:
            raise RuntimeError("Path to genome directory required")


if __name__ == "__main__":
    main()
