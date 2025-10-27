#!/usr/bin/env bash

#$ -pe smp 16
#$ -jc rhel9
#$ -j y
#$ -o drmaa_logs/$JOB_NAME.$JOB_ID
#$ -cwd
#$ -N LXG_phylogeny


if [[ ${CONDA_DEFAULT_ENV} != "comQ_phylogeny" ]]; then
    echo
    echo "This script should be run under the 'comQ_phylogeny' conda environment"
    echo
    echo "The environment can be created with:"
    echo
    echo "     conda env create -f etc/comQ_phylogeny.yaml"
    echo
    echo "then activated with"
    echo
    echo "     conda activate comQ_phylogeny"
    echo
fi

muscle -threads 16 -align NRSW_trees/LXG_tree/for_alignment.fasta -output NRSW_trees/LXG_tree/aligned.fasta

# iqtree run with modelfinder plus , and 1000 bootstrap iterations 
iqtree -s NRSW_trees/LXG_tree/aligned.fasta -m MFP -B 1000 

