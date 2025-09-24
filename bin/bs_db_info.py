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
import gzip
import json
import pickle
import re
import sys

from Bio import SeqIO
from Bio import SearchIO
import click
import pandas as pd
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
        match = matches.group(1)
    else:
        match = 'Unknown'

    return match

def get_spanned_genes(annotation_dir, assembly_accession, seq_accession, start, end):

    """
    Extracts names of genes spanned by hit

    Required params:
        annotation_dir(Path): Path to annotation directory
        assembly_accession(str): GCA accession 
        seq_accession(str): ID of sequence within assembly
        start(str): Start co-ordinate
        end(str): End co-ordinate
    
    Returns:
        genes(str): Comma-separated list of gene names
    """

    gff_file = annotation_dir / assembly_accession / f"{assembly_accession}.gff3"
    gff_cols=('seqid', 'source', 'type', 'start', 'end', 'score', 'strand',
              'phase', 'attributes')
    gff_df = pd.read_csv(gff_file, sep="\t", comment='#', header=0, names=gff_cols, low_memory=False)

    gff = gff_df[(gff_df['seqid'].str.contains(seq_accession)) & (gff_df['type'] == 'CDS') &
              (gff_df['start'] >= int(start)) & (gff_df['end'] <= int(end))]

    gene_re = re.compile(r'gene=([a-zA-Z0-9]*)$')

    genes = gff.apply(lambda x: get_match(x['attributes'], gene_re), axis=1)

    if not genes.empty:
        genes = ','.join([x for x in genes])
        return genes

def parse_blast(genome_path, record, mapping):

    """
    Parses a blast record
    
    Required parameters:
        genome_path(Path): Path to genome directory
        record(pathlib.Path): Path to results file
        mapping(dict): Parsed seq id -> accession info mapping

    Returns:
        parsed_df(pd.DataFrame): DataFrame of parsed results
    """

    parsed_outputs = []
    with record.open() as fh:
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

                genes = get_spanned_genes(genome_path, assembly_accession, seq_accession, start, end)

                description = f"Assembly {assembly_accession}; Strain {strain};"
                if genes:
                    description = f"{description}; Includes genes {genes}"
                hit.description = description

                parsed_results = {
                    'query_id': hit.query_id,
                    'hit_id': hit.id,
                    'pident': pident,
                    'span': span,
                    'mismatches': mismatches,
                    'gaps': gaps,
                    'qstart': qstart,
                    'qend': qend,
                    'hstart': hstart,
                    'hend': hend,
                    'evalue': evalue,
                    'bitscore': bitscore,
                    'description': description
                }
                parsed_outputs.append(parsed_results)

    parsed_df = pd.DataFrame(parsed_outputs)
    return parsed_df

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
@click.option('--output', type=click.Path(exists=False, file_okay=True, dir_okay=False, path_type=Path),
              required=False, help="Path to output file", default=None)
@db.command('tag_info')
def tag_info(db, tags, output):
    """
    Query the database for information on a locus tag
    Reports gene, product, genome accession, organism, sequence id, length,
    coordinates and strand for the assicated gene
    """

    """
    Required arguments:
        db(Path): Path to serialised JSON database
        ids(str): Comma-separated list of locus tags
        output(Path): Path to output file (optional)
    """
    if not Path(db).exists():
        raise FileNotFoundError(f"File '{db}' not found - please build the database first")

    with open(db, 'rb') as fh:
        data = pickle.load(fh)

    if isinstance(tags, str):
        tags = tags.split(',')

    outputs = []
    for locus_id in tags:

        tag = locus_id.split('_')[0]
        if tag in data['locus_tags'].keys():

            locus_data = data['locus_tags'][tag]['genes'][locus_id]

            fields = {
                'Locus tag': locus_id,
                'Genome accession': data['locus_tags'][tag]['accession'],
                'Organism': data['locus_tags'][tag]['organism'],
                'Gene': locus_data['gene'],
                'Product': locus_data['product'],
                'Sequence ID': locus_data['sequence_id'],
                'Sequence length': locus_data['sequence_length'],
                'Coordinates': locus_data['coordinates'],
                'Strand': locus_data['strand']
            }
            outputs.append(fields)

    outputs_df = pd.DataFrame(outputs)

    if len(outputs) == 0:
        raise RuntimeError("No matching records found") 

    if output:
        outputs_df.to_csv(output, sep="\t", index=False)


@click.option('--db', type=click.Path(exists=False, file_okay=True, dir_okay=False, path_type=Path),
              required=False, help="Path to serialised database", default='data/full/bs_genome_info.pkl')
@click.option('--tags', type=str, required=True, help="Comma-separated list of locus tags")
@click.option('--protein_dir', type=click.Path(file_okay=False, dir_okay=True, path_type=Path), required=False,
               help="Path to directory of protein sequences", default="data/full/fasta/proteins")
@click.option('--outfile', type=click.Path(exists=False, file_okay=True, dir_okay=False, path_type=Path),
              required=True, help="Path to output file")

@db.command('fasta')
def fasta(db, tags, protein_dir, outfile):
    """
    Extract protein sequences for a list of locus tags
    """

    """
    Required arguments:
        db(Path): Path to serialised JSON database
        tags(str): Comma-separated list of locus tags
        protein_dir(Path): Path to directory of protein sequences
        outfile(Path): Path to output file
    """

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

@cli.group()
def blast():
    """Commands relating to the blast result handling"""
    pass

@click.option('--genome_path', type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
              required=False, help="Path to genome directory", default='data/full/annotations/')
@click.option('--db', type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=Path),
              required=False, help="Path to serialised database", default='data/full/bs_genome_info.pkl')
@click.option('--blast', type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=Path),
              required=True, help="Path to tab-delimited blast results file")
@click.option('--output', type=click.Path(exists=False, file_okay=True, dir_okay=False, path_type=Path),
              required=True, help="Path to output file", default=None)

@blast.command('genome_hit_info')
def genome_hit_info(genome_path, db, blast, output):
    """
    Parses a blast results file to identify genes spanned by hits, and annotates
    hits with genome accession/strain and any genes spanned by the hit
    """

    """
    Required arguments:
        genome_path(Path): Path to genome directory
        db(Path): Path to serialised JSON database
        blast(Path): Path to tab-delimited blast results file
        output(Path): Path to output file
    """

    if not Path(db).exists():
        raise FileNotFoundError(f"File '{db}' not found - please build the database first")

    with open(db, 'rb') as fh:
        data = pickle.load(fh)

    mapping = data['id_map']

    record = Path(blast)
    parsed_df = parse_blast(genome_path, record, mapping)

    if parsed_df.empty:
        raise RuntimeError("No results found in blast file")

    parsed_df.to_csv(output, sep="\t", index=False)

if __name__ == '__main__':
    cli()