#!/bin/env bash

#$ -cwd
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID
#$ -pe smp 32

gtdbtk classify_wf --genome_dir gtdbtk/genomes --out_dir gtdbtk/results -x fasta  \
                          --cpus 32 --pplacer_cpus 32 --tmpdir $TMPDIR
