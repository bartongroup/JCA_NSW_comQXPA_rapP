#!/bin/env bash

#$ -pe smp 24
#$ -j y
#$ -cwd
#$ -o job_logs/$JOB_NAME.$JOB_ID

mkdir -p results

skani triangle -t 24 --slow -s 50 -o results/skani.out --detailed --full-matrix --diagonal refined/fasta/genomes/*

