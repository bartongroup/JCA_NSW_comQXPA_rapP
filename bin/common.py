'''
Functions reused in multiple scripts live here...
'''

import sys
import re
from json import load

from lxml import etree
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry

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

#def get_busco_excludes(busco_threshold):
#    """
#    Obtains list of accessions which do not meet defined BUSCO completeness threshold
#
#    Required paramsters:
#        busco_threshold (int): Proportion of complete BUSCOs required
#    
#    Returns:
#        exclude(list): list of accessions to exclude
#    """
#
#    busco_results = list(Path('data/full/busco').glob('*/short_summary.specific*json'))
#    exclude = []
#
#    for result in busco_results:
#
#        file = result.name
#        accession = file.replace('short_summary.specific.bacillales_odb10.', '').replace('.json','')
#
#        with open(result, encoding='UTF-8') as fh:
#            dat = load(fh)
#            complete = dat.get('results').get('Complete percentage')
#
#            if float(complete) < busco_threshold:
#                exclude.append(accession)
#
#    return(exclude)

def get_busco_completeness(busco_lineage, accession):

    """
    Obtains busco completeness percentage for a provided genome accession

    Required params: 
        busco_lineage(str): BUSCO lineage used for classification
        accession(str): Genome accession

    Returns:
        completeness(float): BUSCO completeness percentage
    """

    busco_result = f'data/full/busco/{accession}/short_summary.specific.{busco_lineage}.{accession}.json'

    with open(busco_result, encoding='UTF-8') as fh:
        dat = load(fh)
        completeness = dat.get('results').get('Complete percentage')

    return completeness

def join_non_empty(vals):

    """
    Joins strings passed as list with semicolons if both are not empty,
    otherwise returns first or second element if not empty (or second element)
    as fall-back...

    Required params: 
        vals(list): list to join

    Returns
        text(str): joined string, or single value
    """

    if vals[0] and vals[1]:
        return "; ".join(vals)

    if vals[0]:
        return vals[0]

    return vals[1]


def combine_metadata(metadata, classifications, busco_lineage, busco_threshold, output):

    """
    Update metadata sheet to add GTDBTK classifications and BUSCO completeness, 
    and add columns to indicate wheter genome should be retained ('keep') and 
    reason[s] for exclusions

    Required params: 
        metadata(str): Path to metadata file created by 00_build_genome_dbs.py
        classifications(str): Path to GTDBTK classification output file
        busco_lineage(str): BUSCO lineage 
        busco_threshold(int): Proportion of complete BUSCOs required for inclusion
        output(str): Path to output file to create
    
    Returns:
        None
    """
    metadata = pd.read_csv(metadata, sep="\t")
    classifications = pd.read_csv(classifications, sep="\t")

    # Add species column from GTDBTK outputs
    classifications['Species'] = classifications['classification'].map(lambda x: x.split(';')[-1].replace('s__',''))
    #classifications.loc[classifications['Species'] == "Bacillus subtilis", 'Keep'] = 1
    classifications = classifications[['user_genome','Species']]

    metadata = pd.merge(metadata, classifications, how='inner', left_on='Accession', right_on='user_genome')
    metadata = metadata.drop('user_genome', axis = 1)

    # Add BUSCO completeness...
    metadata['BUSCO completeness'] = metadata['Accession'].map(lambda x: get_busco_completeness(busco_lineage, x))

    # Add 'Keep' column defining accessions to be retained, and column explaining rationale for this
    metadata['Keep'] = 1
    metadata['GTDBTK exclusion'] = ""
    metadata['BUSCO exclusion'] = ""

    metadata.loc[metadata['Species'] != "Bacillus subtilis", 'Keep'] = 0
    metadata.loc[metadata['Species'] != "Bacillus subtilis", 'GTDBTK exclusion'] = 'GTDBTK classification'
    metadata.loc[metadata['BUSCO completeness'] < busco_threshold, 'Keep'] = 0
    metadata.loc[metadata['BUSCO completeness'] < busco_threshold, 'BUSCO exclusion'] = 'Below BUSCO threshold'

    metadata["Exclusion criteria"] = metadata[["GTDBTK exclusion", "BUSCO exclusion"]].apply(join_non_empty, axis=1)
    metadata.drop(['GTDBTK exclusion', 'BUSCO exclusion'], axis = 1, inplace = True)

    metadata.to_csv(output, sep = "\t", index = False)
