#!/usr/bin/env python

"""
Tool to extract various information/sequences from the combined database

Use cases include:
    Identifying genes spanned by blast hits from a blast results file
    Reporting gene information/sequence relating to a locus tag

Uses a serialised JSON database of genome information, which will be generated if
not present
"""

from pathlib import Path
import json
import pickle
import re
import sys

from Bio import SeqIO
import click
from tqdm import tqdm

sys.tracebacklimit = 1

def locus_tag_to_GCA(locus_tag, data):
    """
    Converts a locus tag to a GCA format accession  
    Required arguments:
        locus_tag(str): Locus tag to convert
        data(dict): Database of information

    Returns:
        accession(str): Corresponding genome accession
    """

    short_tag = locus_tag.split('_')[0]

    if short_tag in data['locus_tags'].keys():
        accession = data['locus_tags'][short_tag]['genes'][locus_tag]['sequence_id']
        match = re.match(r'lcl\|(GCA_[0-9]+)\|', accession)
        if match:
            accession = match.group(1)
        else:
            raise RuntimeError(f"Could not parse accession from '{accession}'")
    else:
        raise KeyError(f"No accession found for locus tag '{locus_tag}'")

    return accession

def parse_genomes(genome_dir):

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
        genomes = list(genome_dir.glob('*/*embl'))
        compressed = False

    if len(genomes) == 0:
        raise RuntimeError(f"No genomes (.embl/.gz) found in '{genome_dir}'")

    data= {}
    for genome in tqdm(genomes, desc='Parsing genomes...'):

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

        fh.close()

        if short_tag:
            data[short_tag] = genome_data

    return data

@click.group()
def cli():
    """Command line tool to extract information from the combined B.subtilis database"""
    pass

@cli.group()
def db():
    """Commands relating to the serialised database"""
    pass

@click.option('--db', type=click.Path(exists=False, file_okay=True, dir_okay=False, path_type=Path),
              required=False, help="Path to serialised database", default='data/full/bs_genome_info.pkl')
@click.option('--genome_path', type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
              required=False, help="Path to genome directory", default='data/full/annotations/')
@click.option('--id_map', type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=Path),
              required=False, help="Path to id mapping json file", default='data/full/id_mapping.json')
@click.option('--force', is_flag=True, help="Force rebuild of database if already present", default=False)

@db.command('build')
def build(genome_path, id_map, db, force):
    """Build the serialised database from genome records and id mapping file"""

    if force and Path(db).exists():
        Path(db).unlink()

    if Path(db).exists():
        raise FileExistsError(f"File '{db}' already exists - please remove it, or rerun with `--force` if you wish to rebuild")

    data={}
    data['locus_tags'] = parse_genomes(genome_path)

    with open(id_map, 'r', encoding='UTF-8') as fh:
        data['id_map'] = json.load(fh)

    with open(db, 'wb') as out_fh:
        pickle.dump(data, out_fh)

@click.option('--db', type=click.Path(exists=False, file_okay=True, dir_okay=False, path_type=Path),
              required=False, help="Path to serialised database", default='data/full/bs_genome_info.pkl')
@click.option('--tags', type=str, required=True, help="Comma-separated list of locus tags")
@db.command('info')
def info(db, tags):
    """
    Query the database for information on a locus tag

    Required arguments:
        db(Path): Path to serialised JSON database
        ids(str): Comma-separated list of locus tags
    """
    if not Path(db).exists():
        raise FileNotFoundError(f"File '{db}' not found - please build the database first")

    with open(db, 'rb') as fh:
        data = pickle.load(fh)

    if isinstance(tags, str):
        tags = tags.split(',')

    for locus_id in tags:

        tag = locus_id.split('_')[0]
        if tag in data['locus_tags'].keys():

            locus_data = data['locus_tags'][tag]['genes'][locus_id]
            fields = (locus_id, data['locus_tags'][tag]['accession'], data['locus_tags'][tag]['organism'],
                    f"gene: {locus_data['gene']}", locus_data['product'], locus_data['sequence_id'], str(locus_data['sequence_length']),
                    locus_data['coordinates'], str(locus_data['strand']))

            print("\t".join(fields))
        else:
            print(tag)

@click.option('--db', type=click.Path(exists=False, file_okay=True, dir_okay=False, path_type=Path),
              required=False, help="Path to serialised database", default='data/full/bs_genome_info.pkl')
@click.option('--tags', type=str, required=True, help="Comma-separated list of locus tags")
@click.option('--protein_dir', type=click.Path(file_okay=False, dir_okay=True, path_type=Path), required=False,
               help="Path to directory of protein sequences", default="data/full/fasta/proteins")
@click.option('--outfile', type=click.Path(exists=False, file_okay=True, dir_okay=False, path_type=Path),
              required=True, help="Path to output file")

@db.command('fasta')
def fasta(db, tags, protein_dir, outfile):

    with open(db, 'rb') as fh:
        data = pickle.load(fh)

    if isinstance(tags, str):
        tags = tags.split(',')

    output_records=[]
    for tag in tags:
        accession = locus_tag_to_GCA(tag, data)
        with open(Path(protein_dir / f"{accession}.fasta"), 'r', encoding='UTF-8') as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                if tag in record.id:
                    output_records.append(record)
                    break

    if len(output_records) == 0:
        raise RuntimeError("No matching records found")

    with open(outfile, 'w', encoding='UTF-8') as outfh:
        SeqIO.write(output_records, outfh, format='fasta')

if __name__ == '__main__':
    cli()