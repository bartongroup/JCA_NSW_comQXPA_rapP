'''
Functions reused in multiple scripts live here...
'''

def get_taxa_name():

    """
    Looks up taxa via ENA browser API to retrieve scientific name

    Required parameters:
        None

    Returns:
        species_name(str): parsed scientific name
    """

    uri = f"{ENA_URI}browser/api/xml/taxon:{NCBI_TAXID}"
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

def get_gene_sequences(local_genome):

    """
    Parses the downloaded genome record to extract the required sequences

    required parameters:
        local_genome(libpath.Path): Path to local copy of record

    returns:
        gene_seqs(dict): dictionary containing gene and amino acid sequences 
                         as BioSeq objects, keyed on seq ID
    """

    gene_seqs = dict()
    accession = str(local_genome.stem).replace('.embl','')

    with gzip.open(local_genome, 'rt') as fh:
        for record in SeqIO.parse(fh, format = 'embl'):
            organism = record.annotations['organism']

            for feature in record.features:

                # Handle genes by extracting subsequence from record
                if feature.type == 'gene':
                    if 'gene' in feature.qualifiers.keys() and feature.qualifiers.get('gene')[0] in REQUIRED_GENES:

                        gene = feature.extract(record)
                        gene_id = feature.qualifiers.get('gene')[0]

                        id = f'{accession}_{gene_id}_gene'

                        gene.id = id
                        gene.description = f"{organism} {gene_id} gene"

                        gene_seqs[id] = gene

                # CDS feature contain the protein sequence as a 'translation' qualifier
                elif feature.type == 'CDS':

                    if 'gene' in feature.qualifiers.keys() and feature.qualifiers.get('gene')[0] in REQUIRED_GENES:

                        if 'translation' in feature.qualifiers.keys():

                            gene_id = feature.qualifiers.get('gene')[0]
                            prot_seq = feature.qualifiers.get('translation')[0]
                            id = f'{accession}_{gene_id}_protein'

                            protein = SeqRecord(Seq(prot_seq), 
                                                id = id, 
                                                description = f"{organism} {gene_id} protein")

                            gene_seqs[id] = protein

    return(gene_seqs)