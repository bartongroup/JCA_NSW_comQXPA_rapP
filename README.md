# JCA_NSW_comP

PI: Nicola Stanley-Wall

A number of mutations have been identified within *comP* in evolved strains of *Bacillus subtilis*. It is needed to know how stable *comP* is overall within *B. subtilis* genomes, and compare to a related kinase. 

This requires a survey of the available assemblies, consequently retrieving and identifying 'intact' genomes, before downdloading these, converting to fasta, blast indexing each genome, and creating an alias database covering the full set. 

Complete genomes are identified by searching the ENA portal API for assemblies with a provided NCBI taxonomy id. This returns json-formatted results inlucing the SRA accesion of each assembly.

The corresponding SRA XML files are then downloaded, and parsed. Intact genomes are identified on the basis of scaffold number, which is required to be <=5 in order to allow for plasmids. Note that it is not sufficient to look for assemblies with 0 gaps since these are not always correctly reported.

EMBL format gzip compressed records for the assemblies are downloaded into a 'genomes' directory, converted to fasta format ('fasta' directory), and blast indexed ('blast_db' directory). A summary spreadsheet of the downloaded genomes is also written to the project root.

## Building genome databases 

A conda environment (named `genome_dbs`) containing all necessary pre-requisites can be created by running:

`conda env create -f etc/conda.yaml`

This will need to be activated prior to running the script with:

`conda activate genome_dbs`

The `bin/build_genome_dbs.py` script includes an `NCBI_TAXID` constant which defines the taxonomy ID to use for retrieving the genomes. This is set to 1423 (*Bacillus subtilis*)  but can be set to any valid NCBI taxonomy accession. Not that when identifying genomes, all child taxa are included in the search i.e. strains with alternate identifiers.

Once the taxonomy identifier is set, simply run

`bin/build_genome_dbs.py` - it will probably take about half an hour for ~400 genomes.

Note some 404 errors (not found) are likely to be reported during the genome download phase - these will most likely result from assemblies which have metadata available but no completed submission.

### Shell modifications

Searching the resulting blast database will require a larger number of files to be opened simultaneously then permitted by the default shell environment.

Run `ulimit -n 4096` to increase this limit from the default of 1024 prior to running a search. This can also be included in a bash script prior to a blast command.

## Substituting improved genomes

Improved versions of some genomes have been made available, include a small number of additional sequences. This uses the excel spreadsheet produced when building the genome database, and identifies which genomes have improved sequences, removes the old versions and replaces them with the new ones. Any additional sequences are also added at this point.

## Reannotation

To provide consistent annotations, a snakemake workflow enables running bakta annotations on every fasta genome, followed by BUSCO completeness checks. This just requires the `bin/run_workflow.sh` script - the workflow will take care of creating conda environments for each tool and downloading the necessary databases prior to running the analysis. It expects the fasta files to be found in a `fasta` directory as created by `build_genome_dbs.py`, and creates it's outputs in a directies named `annotations` and `busco`.  

The `workflow/Snakefile` sets a `busco_lineage` variable at the top of the file which will require updating if applied to a non-Bacilliacae species.

## Complete genome selection

Genomes which are considered 'complete' based upon a BUSCO completeness score >98% are identified using the `bin/select_complete.py` script. This simply parses the json BUSCO outputs to identify those which have >98% completion, and creates copies in a `complete` subdirectory, including embl-formatted genomes, fasta files and blast databases.

The completeness threshold is defined by the `THRESHOLD` constand at the top of the script, while the NCBI taxonomy identifier used is similarly defined within the NCBI_TAX_ID constant. 

## Sequence extraction

*comP* and *comQ* are still inconsitently annotated, however *comX* is reliably named so used as an anchor to obtain the flanking genes. The `extract_sequences.py` script takes care of this and produced a fasta file for each of the three genes in '`complete/fasta/protein` and blast databases in `complete/blast/protein`. A summary spreadsheet containing info about the genes is also created in `complete/gene_summary.xlsx`