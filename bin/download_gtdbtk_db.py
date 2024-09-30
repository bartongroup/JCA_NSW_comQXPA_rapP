#!/bin/env python

"""
Downloads GTDB-TK database sections in parallel, combines and unpacks
"""

import argparse
from multiprocessing import Pool
from pathlib import Path
import shutil
import subprocess
import tarfile

from bs4 import BeautifulSoup
import requests

from common import make_request

def get_gtdbtk_part_list(release):

    """
    Retrieves data on separate tar sections for gtdbtk database

    Required arguments:
        release(str): GTDB database release number

    Returns:
        uris(list): URIs of tar sections
    """

    uri = f"https://data.gtdb.ecogenomic.org/releases/release{release}/{release}.0/auxillary_files/gtdbtk_package/split_package/"
    soup = BeautifulSoup(make_request(uri), features='lxml')

    links = []
    for link in soup.find_all('a'):
        if link.get('href').startswith('gtdb'):
            part_uri = uri + link.get('href')
            links.append(part_uri)

    return links

def download_file(url):

    """ 
    Downloads tar section using requests 
    
    Required parameters:
        url(str): url of file to download
    
    Returns:
        None
    """

    local_filename = url.split('/')[-1]
    local_filename = Path(f'databases/gtdbtk_download/{local_filename}')

    if not Path(local_filename).exists():
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(local_filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)

def unpack(release):
    """ 
    Combines separate parts into single tar archive and unpacks 

    Required params:
        release (str): release version
    
    Returns:
        None
    """

    parts = sorted(list(Path("databases/gtdbtk_download/").glob(f"gtdbtk_r{release}_data.tar.gz.part*")))

    cmd = ['cat'] + parts

    with open(f"databases/gtdbtk_download/gtdbtk_r{release}_data.tar.gz", 'wb') as fh:
        subprocess.run(cmd, check=True, stdout = fh)

    with tarfile.open(f"databases/gtdbtk_download/gtdbtk_r{release}_data.tar.gz", "r") as handle:
        handle.extractall(path='databases/gtdbtk_download/extract')

    for file in list(Path(f"databases/gtdbtk_download/extract/release{release}").glob("*")):
        shutil.move(str(file), "databases/gtdbtk")

    shutil.rmtree("databases/gtdbtk_download")


def main():
    """ Main process """

    parser = argparse.ArgumentParser(
        prog = "download_gtdbtk_db.py",
        description="Downloads split gtdbtk database, merges and unpacks"
    )
    parser.add_argument('-r', '--release', action='store', dest="release", required=True)

    args = parser.parse_args()

    try:
        Path('databases/gtdbtk').mkdir(exist_ok=True)
    except FileExistsError as e:
        print(e)

    try:
        Path('databases/gtdbtk_download').mkdir(exist_ok=True)
    except FileExistsError as e:
        print(e)

    links = get_gtdbtk_part_list(args.release)

    pool = Pool(6)
    try:
        results = pool.map(download_file, links)
        print(results)
    except requests.exceptions.RequestException as e:
        print(f'Download failed: {e}')

    pool.close()
    pool.join()

    unpack(args.release)


if __name__ == "__main__":
    main()
