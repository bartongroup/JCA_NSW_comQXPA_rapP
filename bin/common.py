'''
Functions reused in multiple scripts live here...
'''
from tqdm import tqdm
from Bio import SeqIO
from lxml import etree
from pathlib import Path
import requests
import subprocess
import sys
import re

NCBI_TAXID = '1423'

"""
Identifies all complete B.subtilis genomes available via ENA and downloads them
"""

def make_request(uri):

    """
    Makes HTTP request and returns result body

    Required params:
        uri(str): URI for request

    Returns:
        text of response
    """
    try:
        r = requests.get(uri)
    except Exception as e:
        print(f'Error requesting {uri}: {e}')

    return(r.text)


def get_taxa_name(ncbi_taxid):

    """
    Looks up taxa via ENA browser API to retrieve scientific name

    Required parameters:
        ncbi_taxid(str): Taxonomy identifier of species

    Returns:
        species_name(str): parsed scientific name
    """

    uri = f"https://www.ebi.ac.uk/ena/browser/api/xml/taxon:{ncbi_taxid}"
    result = make_request(uri)

    # lxml does not support embedded encodings...
    result = result.replace(' encoding="UTF-8"', '')
    root = etree.fromstring(result)

    error = root.xpath("/ErrorDetails/status")
    if error:
        print(f"Error obtaining taxonomy lookup")
        print(root)
        sys.exit(1)

    taxon = root.xpath("/TAXON_SET/taxon")[0]
    species_name = taxon.get('scientificName')

    return(species_name)

def blast_index(fasta_dir, blast_dir, species):

    """
    Creates a blast index for each fasta file contained within directory

    Required parameters:
        fasta_dir(Path): Location of fasta files
        blast_dir(Path): Location of blast database

    Returns:
        None
    """

    accessions = list()
    db_stats = {'count': 0, 'length': 0}
    
    fasta_files = fasta_dir.glob('*.fasta')
    for fasta_file in tqdm(fasta_files, desc = 'Blast index'):
        with open(fasta_file, 'r') as fh:
            records = SeqIO.parse(fh, format = 'fasta')
            for record in records:
                db_stats['count'] += 1
                db_stats['length'] += len(record.seq)

        accession = fasta_file.name.replace('.fasta','')
        accessions.append(accession)

        cmd = ["makeblastdb",  "-in",  fasta_file, "-dbtype", "nucl", "-out", blast_dir / Path(accession),
               "-taxid", NCBI_TAXID, "-title", accession]
        try:
            result = subprocess.run(cmd, check = True, capture_output = True)
        except subprocess.CalledProcessError as e:
            print(f"Index failed: {accession}")
            print(result.stdout)
            print(result.stderr)
    
    accessions = ' '.join([f"{a}" for a in accessions])

    index = f'''\
TITLE {species} complete genomes
DBLIST {accessions}
NSEQ {db_stats['count']}
LENGTH {db_stats['length']}
'''.strip()

    with open(f'{blast_dir}/{species.replace(' ','_')}_complete_genomes.nal', 'w') as fh:
        fh.writelines(index)

def clean_description(description):

    """
    tidies up variously formatted descriptions into something reasonably consistent

    Required params:
        description(str): record description
    
    Returns:
        description(str): cleaned description
    """

    description = re.sub(r'[, ]*complete genome[.]*', '', description)
    description = re.sub(r'[, ]*genome assembly[,:_ A-Za-z0-9]*', '', description)
    description = re.sub(r'chromosome[ ,]*', '', description)
    description = re.sub('whole genome shotgun sequence.','', description)
    description = re.sub(r'NODE[A-Za-z0-9_,]*', '', description)
    description = re.sub(r'[ ,]*$','', description)

    return(description)

def clean_strain(strain):

    """
    tidies up variously formatted strain names into something reasonably consistent

    Required params:
        strain(str): strain name
    
    Returns:
        strain(str): cleaned strain name
    """

    strain = re.sub(r'[, ]*(Mutant )?Bacillus subtilis (strain|isolate) [.]*', '', strain)
    strain = re.sub(r'[, ]*(Mutant )?Bacillus subtilis [.]*', '', strain)
    strain = re.sub(r' (strain|isolate)', '', strain)
    #strain = strain.replace(' ', '_')

    return(strain)