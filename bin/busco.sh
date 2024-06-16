#!/bin/env bash

USAGE=$(cat << END
Annotates fasta formatted genome with bakta
Usage: $0 -g genome_name -l lineage -o output_directory"
END
)

while getopts "g:l:o:h" opt; do
    case ${opt} in
        g)
            GENOME=$OPTARG
            ;;
        l)
            LINEAGE=$OPTARG
            ;;
        o) 
            OUTDIR=$OPTARG
            ;;
        h)
            echo $USAGE
            ;;
        ?)
            echo "Invalid option: ${opt}"
            exit 1
            ;;
    esac
done

if [[ -z "${GENOME}" ]] | [[ -z "${LINEAGE}" ]] | [[ -z "${OUTDIR}" ]]; then
    echo ${USAGE}
    exit 1
fi

busco -f -l ${LINEAGE} --mode genome -i fasta/${GENOME}.fasta -c 8 -o ${OUTDIR}




