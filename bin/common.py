'''
Functions reused in multiple scripts live here...
'''

from pathlib import Path
import subprocess
import sys
import re

from Bio import SeqIO
from lxml import etree
import requests
from requests.adapters import HTTPAdapter, Retry
from tqdm import tqdm

def make_request(uri):

    """
    Makes HTTP request and returns result body

    Required params:
        uri(str): URI for request

    Returns:
        text of response
    """

    s = requests.Session()

    retries = Retry(total=10,
              backoff_factor=2,
              status_forcelist=[ 500, 502, 503, 504 ],
              raise_on_status=True)

    s.mount('https://', HTTPAdapter(max_retries=retries))

    try:
        r = s.get(uri, timeout=30)
    except requests.exceptions.RequestException as errex:
        print(f"Exception: {uri} - {errex}")

    return r.text


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
        print("Error obtaining taxonomy lookup")
        print(root)
        sys.exit(1)

    taxon = root.xpath("/TAXON_SET/taxon")[0]
    species_name = taxon.get('scientificName')

    return species_name

def blast_index(fasta_dir, blast_dir, species, ncbi_taxid, index_type):

    """
    Creates a blast index for each genome fasta file contained within directory

    Required parameters:
        fasta_dir(Path): Location of fasta files
        blast_dir(Path): Location of blast database
        species(str): Species name
        ncbi_taxid(str): NCBI taxonomy ID
        index_type(str) : index type (nucl/prot)

    Returns:
        None
    """

    species_name=species.replace(' ','_')

    if index_type not in ['nucl','prot']:
        raise ValueError(f'Invalid index type provided ({index_type})')

    accessions = []
    db_stats = {'count': 0, 'length': 0}

    if index_type == 'nucl':
        blast_dir = blast_dir / 'genomes'
        title = 'TITLE {species} complete genomes'
        meta_db = f'{blast_dir}/{species_name}_complete_genomes.nal'
    else:
        blast_dir = blast_dir / 'proteins'
        title  =  f"TITLE {species} complete proteome"
        meta_db = f'{blast_dir}/{species_name}_proteins.pal'

    fasta_files = fasta_dir.glob('*.fasta')
    for fasta_file in tqdm(fasta_files, desc = 'Blast index'):
        with open(fasta_file, 'r', encoding='UTF-8') as fh:
            records = SeqIO.parse(fh, format = 'fasta')
            for record in records:
                db_stats['count'] += 1
                db_stats['length'] += len(record.seq)

        accession = fasta_file.name.replace('.fasta','')
        accessions.append(accession)

        cmd = ["makeblastdb",  "-in",  fasta_file, "-dbtype", index_type, "-out", blast_dir / Path(accession),
               "-taxid", str(ncbi_taxid), "-title", accession]
        try:
            res = subprocess.run(cmd, check = True, capture_output = False)
        except subprocess.CalledProcessError as e:
            print(f"Index failed: {accession} - {e}")

    accessions = ' '.join([f"{a}" for a in accessions])

    index = f'''\
{title} 
DBLIST {accessions}
NSEQ {db_stats['count']}
LENGTH {db_stats['length']}
'''.strip()

    with open(meta_db, 'w', encoding='UTF-8') as fh:
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

    return description

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

    return strain
