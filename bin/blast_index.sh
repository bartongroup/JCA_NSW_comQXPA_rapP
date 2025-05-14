#!/usr/bin/env bash

set -e

# Generates blast index for each proteome, including overarching alias database

usage() { echo "$0 -f fasta_dir -d blast_db_dir -s species -i taxonomy_id -t [n|p]" 2>&1;  exit 1; }

while getopts "f:d:s:i:t:" opt; do
    case "${opt}" in
        f)
            fasta_dir=${OPTARG}
            ;;
        d)
            blast_dir=${OPTARG}
            ;;
        s)
            species=${OPTARG}
            species_name=$(echo ${OPTARG} | sed s'/ /_/g')
            ;;
        i)
            taxid=${OPTARG}
            ;;
        t)
            if [ "${OPTARG}" == "p" ]; then
                type="prot"
            elif [ "${OPTARG}" == "n" ]; then
                type="nucl"
            else
                echo "Type must be [p]rotein or [n]ucleotide"
                usage
            fi
            ;;
        *)
            usage
            ;;
    esac
done

if [ -z "${fasta_dir}" ] || [ -z "${blast_dir}" ] || [ -z "${species}" ] || \
   [ -z "${taxid}" ] || [ -z "${type}" ]; then
    usage
fi

db_seqs=0
db_length=0
accessions=""

if [ "${type}" == "nucl" ]; then
    blast_dir="${blast_dir}/genomes"
    meta_db="${blast_dir}/${species_name}_complete_genomes.nal"
    title="TITLE ${species} complete genomes"
else
    blast_dir="${blast_dir}/proteins"
    meta_db="${blast_dir}/${species_name}_proteins.pal"
    title="TITLE ${species} complete proteome"
fi

for file in $(ls ${fasta_dir}/*.fasta) ; do

    seqs=$(grep -c '>' $file)
    length=$(grep -v '>' $file|wc -c)

    db_seqs=$((${db_seqs}+${seqs}))
    db_length=$((${db_length}+${length}))

    accession=$(basename $file|sed 's/.fasta//')
    accessions+="${accession} "

    makeblastdb -in $file -dbtype ${type} -out ${blast_dir}/${accession} -taxid ${taxid} -title ${accession}
done

accession_str="${accessions[*]}"

# Generate nal/pal database
cat  << EOF > $meta_db
${title}
DBLIST ${accession_str}
NSEQ ${db_seqs}
LENGTH ${db_length}
EOF