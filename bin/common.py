'''
Functions reused in multiple scripts live here...
'''

from pathlib import Path
import subprocess
import sys
import re

from Bio import SeqIO
from lxml import etree
import pandas as pd
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

def combine_metadata(metadata, classifications, output):

    """
    Update metadata sheet to add GTDBTK classifications and BUSCO completeness, 
    and add columns to indicate wheter genome should be retained ('keep') and 
    reason[s] for exclusions

    Required params: 
        metadata(str): Path to metadata file created by 00_build_genome_dbs.py
        classifications(str): Path to GTDBTK classification output file
        output(str): Path to output file to create
    
    Returns:
        None
    """
    metadata = pd.read_csv(metadata, sep="\t")
    classifications = pd.read_csv(classifications, sep="\t")

    # Add species column, and 'keep' column, set to 1 where species is subtilis
    classifications['Species'] = classifications['classification'].map(lambda x: x.split(';')[-1].replace('s__',''))
    classifications['Keep'] = 0
    classifications['Exclusion criteria'] = ""
    classifications.loc[classifications['Species' ]== "Bacillus subtilis", 'Keep'] = 1
    classifications.loc[classifications['Species' ]!= "Bacillus subtilis", 'Exclusion criteria'] = 'GTDBTK classification'
    classifications = classifications[['user_genome','Species', 'Keep', 'Exclusion criteria']]

    metadata = pd.merge(metadata, classifications, how='inner', left_on='Accession', right_on='user_genome')
    metadata = metadata.drop('user_genome', axis = 1)
    metadata.to_csv(output, sep="\t", index=False)