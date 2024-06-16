#!/bin/env bash

USAGE=$(cat << END
Annotates fasta formatted genome with bakta
Usage: $0 -n genome_name -g genus -s species -o output_directory"
END
)

while getopts "n:g:s:o:h" opt; do
    case ${opt} in
        n)
            GENOME=$OPTARG
            ;;
        s)
            SPECIES=$OPTARG
            ;;
        g)
            GENUS=$OPTARG
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

if [[ -z "${GENOME}" ]] | [[ -z "${GENUS}" ]] | [[ -z "${SPECIES}" ]] | [[ -z "${OUTDIR}" ]]; then
    echo ${USAGE}
    exit 1
fi

cp -v fasta/${GENOME}.fasta ${TMPDIR}

bakta --db db --prefix ${GENOME} --genus ${GENUS} --species ${SPECIES} --compliant \
  --threads 8 --output ${TMPDIR}/bakta ${TMPDIR}/${GENOME}.fasta

cp -v ${TMPDIR}/bakta/* ${OUTDIR}