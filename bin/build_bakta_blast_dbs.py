#!/bin/env python

"""
Creates per-genome blast databases from bakta outputs for cds and proteins
"""

from pathlib import Path
from shutil import copy
import subprocess
from tqdm import tqdm
from Bio import SeqIO
from common import get_taxa_name

NCBI_TAXID = "1423"

def main():

    species = get_taxa_name(NCBI_TAXID)

    cds_fastas = Path('complete' , 'fasta', 'cds')
    pep_fastas = Path('complete' , 'fasta', 'protein')

    cds_blasts = Path('complete' , 'blast', 'cds')
    pep_blasts = Path('complete' , 'blast', 'protein')


    cds_fastas.mkdir(exist_ok = True)
    pep_fastas.mkdir(exist_ok = True)
    cds_blasts.mkdir(exist_ok = True)
    pep_blasts.mkdir(exist_ok = True)

    genomes = list(Path('annotations').glob('*'))

    accessions = list()
    db_stats = {'cds_count': 0, 'cds_length': 0,
                'pep_count': 0, 'pep_length': 0}

    for genome in tqdm(genomes):
        accession = genome.name
        accessions.append(accession)

        cds_fasta = genome / f"{accession}.ffn"
        pep_fasta = genome / f"{accession}.faa"

        out_cds_fasta = cds_fastas / f"{accession}.ffn"
        out_pep_fasta = pep_fastas / f"{accession}.ffn"
        cds_blast = cds_blasts / accession
        pep_blast = pep_blasts / accession

        # rewrite ids to include genome accession
        out_records = list()
        with open(cds_fasta, 'r') as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                record.id = f"{accession}-{record.id}"
                out_records.append(record)
                db_stats['cds_count'] += 1
                db_stats['cds_length'] += len(record.seq)
        
        with open(out_cds_fasta, 'w') as fh:
            SeqIO.write(out_records, fh, 'fasta')

        out_records = list()
        with open(pep_fasta, 'r') as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                record.id = f"{accession}:{record.id}"
                out_records.append(record)
                db_stats['pep_count'] += 1
                db_stats['pep_length'] += len(record.seq)
        
        with open(out_pep_fasta, 'w') as fh:
            SeqIO.write(out_records, fh, 'fasta')

        cmd = ["makeblastdb",  "-in",  cds_fasta, "-dbtype", "nucl", "-out", cds_blast,
            "-taxid", NCBI_TAXID, "-title", accession]
        try:
            result = subprocess.run(cmd, check = True, capture_output = True)
        except subprocess.CalledProcessError as e:
            print(f"Index failed: {accession}")
            print(result.stdout)
            print(result.stderr)

        cmd = ["makeblastdb",  "-in", pep_fasta, "-dbtype", "prot", "-out", pep_blast,
            "-taxid", NCBI_TAXID, "-title", accession]
        try:
            result = subprocess.run(cmd, check = True, capture_output = True)
        except subprocess.CalledProcessError as e:
            print(f"Index failed: {accession}")
            print(result.stdout)
            print(result.stderr)

    accessions = ' '.join([f"{a}" for a in accessions])

    index = f'''\
TITLE {species} CDS sequences
DBLIST {accessions}
NSEQ {db_stats['cds_count']}
LENGTH {db_stats['cds_length']}
'''.strip()

    with open(f'{cds_blasts}/{species.replace(' ','_')}_cds.nal', 'w') as fh:
        fh.writelines(index)
    
    index = f'''\
TITLE {species} protein sequences
DBLIST {accessions}
NSEQ {db_stats['pep_count']}
LENGTH {db_stats['pep_length']}
'''.strip()

    with open(f'{pep_blasts}/{species.replace(' ','_')}_proteins.pal', 'w') as fh:
        fh.writelines(index)


if __name__ == "__main__":
    main()