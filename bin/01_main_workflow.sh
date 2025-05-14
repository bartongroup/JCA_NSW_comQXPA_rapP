#!/bin/env bash

# Runs a Snakemake workflow to annotate downloaded genomes with bakta,
# assess completeness with busco and classify taxonomy with GTDB-TK

set -e 

# Some sanity checking...
# work out where we are running from, taking symlinks into account
run_dir=$(readlink -f "$0")
analysis_root=$(dirname ${run_dir}|sed 's/\/bin//')

# now work out where we are, taking symlinks into account
#'curr_dir=$(pwd|sed 's/\/analysis\/[A-Za-z0-9]*$//')
curr_dir=$(pwd)
real_base_dir=$(readlink -f $curr_dir)
if [[ ! -z "${real_base_dir}" ]]; then
	base_dir=${real_base_dir}
fi

# Are we in the right place?
if [[ "${base_dir}" != "${analysis_root}" ]]; then
	echo "This script should be run from the top-level directory of the project"
	exit 1
fi

snakemake --profile profile $@

#snakemake  --snakefile workflow/Snakefile --rerun-incomplete --cluster-config config/workflow.json \
#    --latency-wait 60 -j 50 --jn annotate_{jobid} --use-conda --conda-frontend mamba -k --drmaa \
#	" -cwd -j y -o logs/{rule}.{wildcards}.log -jc {cluster.jc} -pe smp {cluster.threads} -mods l_hard mfree {cluster.mem_free} -adds l_hard local_free {cluster.local_free} " $@