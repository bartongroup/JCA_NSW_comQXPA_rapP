#!/usr/bin/env bash

#$ -cwd
#$ -o job_logs/$JOB_NAME.$JOB_ID
#$ -j y
#$ -pe smp 48 
#$ -jc 4week

set -e 

cp -RLv final/fasta/genomes $TMPDIR
# add the outgroup
cp -v refined/fasta/genomes/GCA_000772125.fasta $TMPDIR/genomes

RAXML=$(which raxmlHPC-PTHREADS)

GENOMES=$(ls ${TMPDIR}/genomes/*.fasta)

GENOME_ARGS=""

for GENOME in ${GENOMES[@]}; do
    ACCESSION=$(basename $GENOME|sed 's/.fasta//')
    GENOME_ARGS+="--genome ${ACCESSION}=${GENOME} "
done

bcgTree.pl --threads 48 --raxml-bin=${RAXML} --raxml-aa-substitution-model GTR --outdir $TMPDIR/bcgTree ${GENOME_ARGS} 
cp -Rv $TMPDIR/bcgTree results