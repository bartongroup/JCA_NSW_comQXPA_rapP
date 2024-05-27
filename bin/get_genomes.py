#!/bin/env python

import requests
from json import loads

"""
Identifies all complete B.subtilis genomes available via ENA and downloads them
"""

SEARCH_URI="https://www.ebi.ac.uk/ena/portal/api/search"

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

def search_available():
    """
    Carries out search of ENA using portal api to identify available genomes 

    required parameters: 
        None

    returns:
        genome_info(dict)
    """

    uri = f"{SEARCH_URI}/?result=assembly&query=tax_tree(1423)&fields=assembly_title%2Ctax_id&format=json"
    result = make_request(uri)
    json = loads(result)

    return(result)


def main():

    genome_info = search_available()
    print(genome_info)


if __name__ == "__main__":
    main()