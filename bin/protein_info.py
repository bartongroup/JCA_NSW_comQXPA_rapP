#!/usr/bin/env python

from argparse import ArgumentParser
from pathlib import Path
import gzip
import json
import re
import sys
from io import TextIOWrapper

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

    genomes = list(genome_dir.glob('*/*embl'))

    if len(genomes) == 0:
        raise RuntimeError(f"No genomes (.embl) found in '{genome_dir}'")

    data= {}
    for genome in genomes:

        genome_data = {}
        genome_data['genes'] = {}
        locus_tag = ""
        short_tag = ""

        with open(genome, encoding = 'UTF-8') as fh:

            for record in SeqIO.parse(fh, format = 'embl'):
                genome_accession = genome.name
                sequence_id = record.id
                sequence_length = len(record.seq)
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
                        genome_data['genes'][locus_tag]['sequence_length'] = sequence_length
                        genome_data['genes'][locus_tag]['coordinates'] = f"{feature.location.start}-{feature.location.end}"
                        genome_data['genes'][locus_tag]['strand'] = feature.location.strand

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

    if isinstance(ids, str):
        ids= ids.split()

    for locus_id in ids:

        tag = locus_id.split('_')[0]
        if tag in data:

            locus_data = data[tag]['genes'][locus_id]
            fields = (locus_id, data[tag]['accession'], data[tag]['organism'],
                    f"gene: {locus_data['gene']}", locus_data['product'], locus_data['sequence_id'], str(locus_data['sequence_length']),
                    locus_data['coordinates'], str(locus_data['strand']))

            print("\t".join(fields))
        else:
            print(tag)

def print_fasta(ids, data, protein_dir, outfile):

    """
    Looks up and prints fasta-formatted sequence[s] for
    provided identifiers

    Required params:
        ids(string): (posibly multiline) list of ids
        data(dict): Dictionary keyed on locus tax prefix produced by parse_data()
        protein_dir(pathlib.Path): Path to protein fasta directory
        outfile (str): Path to output file
    
    Returns:
        None
    """

    with open(outfile, 'w', encoding='UTF-8') as outfh:
        if isinstance(ids, str):
            ids= ids.split()

        for locus_id in ids:

            tag = locus_id.split('_')[0]
            if tag in data:
                accession = data[tag]['accession']

                with open(Path(protein_dir / f"{accession}.fasta"), 'r', encoding='UTF-8') as fh:
                    for record in SeqIO.parse(fh, 'fasta'):
                        if locus_id in record.id:
                            SeqIO.write([record], outfh, format='fasta')
                            break


def main():

    """
    Main method...
    """

    parser = ArgumentParser(
        prog = __file__,
            description = """
            Looks up genome/gene associated with a locus tagged identifier i.e.
            HLBEOP_22245, or returns the appropriate fasta formatted sequence. 
            If preparsed data is available (locus_tag_data.json)
            this will be used, otherwise the available genomes will be parsed, 
            and the results stored for reuse""" )

    parser.add_argument('-i', '--input', help="Idenfiers(s) to resolve")
    parser.add_argument('-f', '--fasta', required=False, action='store_true', help="Return fasta formatted sequence")
    parser.add_argument('-o', '--output', required=False, help='Path to write output fasta file')
    parser.add_argument('-j', '--json', help="Path to genome_locus_tags.json.gz file",
        default = 'data/full/genome_locus_tags.json.gz')
    args = parser.parse_args()

    serialised_data = Path(args.json)
    genome_dir = serialised_data.parent / 'annotations'
    protein_dir = serialised_data.parent / 'fasta/proteins/'

    if serialised_data.exists():

        with gzip.open(serialised_data, 'rb') as fh:
            json_data = fh.read()
        data = json.loads(json_data)

    else:
        data = parse_data(genome_dir, serialised_data)

    if args.input:
        with open(args.input, 'r', encoding = 'UTF-8') as fh:
            ids = fh.readlines()
            ids = [id.strip() for id in ids]

    else:
        ids = sys.stdin.read()

    if args.fasta:
        print_fasta(ids, data, protein_dir, args.output)
    else:
        print_details(ids, data)

if __name__ == "__main__":
    main()