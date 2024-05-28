# JCA_NSW_comP

PI: Nicola Stanley-Wall

A number of mutations have been identified within *comP* in evolved strains of *Bacillus subtilis*. It is needed to know how stable *comP* is overall within *B. subtilis* genomes, and compare to a related kinase. 

This requires a survey of the available assemblies to download intact builds and extract the required genes

The `bin/get_genes.py` script defines a `REQUIRED_GENES` list constant at the top of the script which should include the gene names to be retrieved. The NCBI_TAXID constant also readily allows the script to be repurposed for other organsimsm and genes.

An EBI search is initially carried out to identify all assemblies, and the relevent SRA assembly XML document retrieved. This is then parsed and the assembly attributes checked to identify sequences with 0 spanned and unspanned gaps i.e. intact assemblies.

Genome sequences are downloaded and stored as gzipped files in the `genomes` directory, and subsequently parsed to extract the genes and proteins of interest which are written to separate files in the `outputs` directory