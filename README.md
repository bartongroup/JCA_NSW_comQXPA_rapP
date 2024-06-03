# JCA_NSW_comP

PI: Nicola Stanley-Wall

A number of mutations have been identified within *comP* in evolved strains of *Bacillus subtilis*. It is needed to know how stable *comP* is overall within *B. subtilis* genomes, and compare to a related kinase. 

This requires a survey of the available assemblies, consequently retrieving and identifying 'intact' genomes, before downdloading these, converting to fasta, blast indexing each genome, and creating an alias database covering the full set. 

Complete genomes are identified by searching the ENA portal API for assemblies with a provided NCBI taxonomy id. This returns json-formatted results inlucing the SRA accesion of each assembly.

The corresponding SRA XML files are then downloaded, and parsed. Intact genomes are identified on the basis of scaffold number, which is required to be <=5 in order to allow for plasmids. Note that it is not sufficient to look for assemblies with 0 gaps since these are not always correctly reported.

EMBL format gzip compressed records for the assemblies are downloaded into a 'genomes' directory, converted to fasta format ('fasta' directory), and blast indexed ('blast_db' directory). A summary spreadsheet of the downloaded genomes is also written to the project root.

## Usage

A conda environment (named `genome_dbs`) containing all necessary pre-requisites can be created by running:

`conda env create -f etc/conda.yaml`

This will need to be activated prior to running the script with:

`conda activate genome_dbs`

The `bin/build_genome_dbs.py` script includes an `NCBI_TAXID` constant which defines the taxonomy ID to use for retrieving the genomes. This is set to 1423 (*Bacillus subtilis*)  but can be set to any valid NCBI taxonomy accession. Not that when identifying genomes, all child taxa are included in the search i.e. strains with alternate identifiers.

Once the taxonomy identifier is set, simply run

`bin/build_genome_dbs.py` - it will probably take about half an hour for ~400 genomes.

Note some 404 errors (not found) are likely to be reported during the genome download phase - these will most likely result from assemblies which have metadata available but no completed submission.