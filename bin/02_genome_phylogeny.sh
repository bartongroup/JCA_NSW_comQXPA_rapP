#!/usr/bin/env bash

#$ -pe smp 16
#$ -jc rhel9
#$ -j y
#$ -o drmaa_logs/$JOB_NAME.$JOB_ID
#$ -cwd

# Generates genome phylogeny using initial stages of bcgTree,
# but uses IQtree in place of RaxML for performance reasons.
# This is run as a separate script in it's own conda environment 
# to allow the bcgTree.pl script to be modified to skip the RaxML 
# tree building

if [[ ${CONDA_DEFAULT_ENV} != "bsubtilis_phylogeny" ]]; then
    echo
    echo "This script should be run under the 'bsubtilis_phylogeny' conda environment"
    echo
    echo "The environment can be created with:"
    echo
    echo "     conda env create -f etc/bsubtilis_phylogeny.yaml"
    echo
    echo "then activated with"
    echo
    echo "     conda activate bsubtilis_phylogeny"
    echo
fi

# Comment out line which calls raxml in bcgTree.pl
sed -i.bak 's/$bcgTree->run_raxml/#$bcgTree->run_raxml/' $CONDA_PREFIX/bin/bcgTree.pl

# Create bcgTree concatenated alignment from proteomes selected 
# by refine_genome_selection snakemake rule
bcgTree.pl @results/phylogeny/bcgtree.config --outdir results/phylogeny

# iqtree run using an edge proportional parition model (-p) with modelfinder
# plus allowing partition merging, and 1000 bootstrap iterations using
# resampling of partitions followed by sites within partitions 
# (--sampling GENESITE)

iqtree -s results/phylogeny/full_alignment.concat.fa -p results/phylogeny/full_alignment.concat.partition -m MFP+MERGE -B 1000 --sampling GENESITE -nt AUTO -ntmax 16

