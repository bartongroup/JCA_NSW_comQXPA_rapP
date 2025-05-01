#!/usr/bin/env bash

#$ -cwd
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID
#$ -pe smp 32

GENOME_DIR="refined/fasta/genomes"
mkdir -p results/fastANI

ls $GENOME_DIR/*.fasta > results/fastANI/genomes.txt

fastANI --ql results/fastANI/genomes.txt --rl results/fastANI/genomes.txt -t 32 -o results/fastANI/output.txt