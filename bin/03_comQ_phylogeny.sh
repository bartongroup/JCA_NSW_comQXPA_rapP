#!/usr/bin/env bash

#$ -pe smp 16
#$ -jc rhel9
#$ -j y
#$ -o drmaa_logs/$JOB_NAME.$JOB_ID
#$ -cwd
#$ -N comQ_phylogeny

# Generates comQ phylogeny using from refined sequence selection

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

cat data/refined/fasta/target_proteins/comQ/* > $TMPDIR/query.fasta
cat reference_sequences/outgroup/target_proteins/comQ.fasta >> $TMPDIR/query.fasta
sed -i 's/_comQ//g' $TMPDIR/query.fasta
mkdir -p results/comQ_phylogeny

muscle -threads 16 -align $TMPDIR/query.fasta -output results/comQ_phylogeny/comQ.fasta

# iqtree run with modelfinder plus , and 1000 bootstrap iterations 
iqtree -s results/comQ_phylogeny/comQ.fasta -m MFP -B 1000 

