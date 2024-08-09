#!/bin/env bash

#$ -cwd
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID
#$ -pe smp 32

mkdir -p results/gtdbtk
cd results/gtdbtk
ln -s ../../refined/fasta/genomes
cd -

gtdbtk classify_wf --genome_dir results/gtdbtk/genomes --out_dir results/gtdbtk/b_subtilis -x fasta  \
                   --skip_ani_screen --cpus 32 --pplacer_cpus 32 \
                   --tmpdir $TMPDIR --force
