#!/usr/bin/env python

from pathlib import Path
from Bio import SeqIO

genome_dir = Path('final', 'genomes')
genomes = genome_dir.glob('*.embl')
accessions = [str(x.name).replace('.embl', '') for x in genomes]

# Add the outgroup....
accessions.append('GCA_000772125')

out_records = list()

with open('operon/operons.fasta') as fh:
    for record in SeqIO.parse(fh, 'fasta'):
        if record.id in accessions:
            out_records.append(record)
        

with open('operon/operons_final.fasta', 'w') as fh:
    SeqIO.write(out_records, fh, 'fasta')