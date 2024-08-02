#!/bin/env bash

#$ -cwd
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID
#$ -pe smp 32

GENOME_DIR="refined/fasta/genomes"
mkdir -p refined/fastANI

ls $GENOME_DIR/*.fasta > refined/fastANI/genomes.txt

fastANI --ql refined/fastANI/genomes.txt --rl refined/fastANI/genomes.txt -t 32 -o refined/fastANI/output.txt