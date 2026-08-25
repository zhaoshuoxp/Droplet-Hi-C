#!/usr/bin/env bash

set -Eeuo pipefail

threads=24
res_array=(10)
sample=""

while getopts "s:r:t:" opt; do
  case ${opt} in
    s) sample=${OPTARG} ;;
    r) IFS=',' read -ra res_array <<< "${OPTARG}" ;;
    t) threads=${OPTARG} ;;
    *) echo "Usage: $0 [-s sample] [-r res1,res2,...] [-t threads]"; exit 1 ;;
  esac
done

if [ -z "${sample}" ]; then
  echo "Sample name (-s) is required."
  exit 1
fi

if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Threads must be a positive integer."
    exit 1
fi

for r in "${res_array[@]}"; do
    if [[ ! "${r}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Resolution values must be positive integers."
        exit 1
    fi
done

current=$(pwd)
imputed_matrix="${current}/hicluster/imputed_matrix"

for r in "${res_array[@]}"; do
    echo "Processing resolution ${r}kb..."

    cell_table="${imputed_matrix}/${r}kb_resolution/filelist/cell_table.tsv"

    if [[ ! -s "${cell_table}" ]]; then
        echo "Cell table not found for ${r}kb at ${cell_table}" >&2
        exit 1
    fi

    mkdir -p "${imputed_matrix}/${r}kb_resolution/domain/"

    hicluster domain \
        --cell_table_path "${cell_table}" \
        --output_prefix "${imputed_matrix}/${r}kb_resolution/domain/${sample}" \
        --resolution "${r}000" \
        --window_size 20 \
        --cpu "${threads}"

    echo "Domain calling done for ${r}kb"
done
