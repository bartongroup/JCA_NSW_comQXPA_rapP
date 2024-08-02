#!/bin/env bash

#$ -pe smp 24
#$ -j y
#$ -cwd
#$ -o job_logs/$JOB_NAME.$JOB_ID

skani triangle -t 24 -s 50 -o refined/skani.out --detailed --full-matrix --diagonal refined/fasta/genomes/*

