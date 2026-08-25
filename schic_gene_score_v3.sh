#!/usr/bin/env bash

set -Eeuo pipefail

genome="mm10"
threads=24
res_array=(10)
sample=""

while getopts "s:g:r:t:" opt; do
  case ${opt} in
    s) sample=${OPTARG} ;;
    g) genome=${OPTARG} ;;
    r) IFS=',' read -ra res_array <<< "${OPTARG}" ;;
    t) threads=${OPTARG} ;;
    *) echo "Usage: $0 [-s sample] [-g genome] [-r res1,res2,...] [-t threads]"; exit 1 ;;
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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOME_DIR="${DROPLET_HIC_GENOME_DIR:-/home/quanyiz/genome}"

if [[ "${genome}" != "mm10" && "${genome}" != "hg38" ]]; then
    echo "Unsupported genome: ${genome}"
    exit 1
fi

gene_meta="${SCRIPT_DIR}/01.pre-process/bed/${genome}_gene_anno.bed"
chrom_size="${GENOME_DIR}/${genome}/${genome}_noYM.chrom.sizes"

for required_file in "${gene_meta}" "${chrom_size}"; do
    if [[ ! -r "${required_file}" ]]; then
        echo "Required input is not readable: ${required_file}" >&2
        exit 1
    fi
done

current=$(pwd)
imputed_matrix="${current}/hicluster/imputed_matrix"

for r in "${res_array[@]}"; do
    echo "Processing resolution ${r}kb..."

    mkdir -p "${imputed_matrix}/${r}kb_resolution/filelist"
    cell_table="${imputed_matrix}/${r}kb_resolution/filelist/cell_table.tsv"

    echo "Generating/Updating cell table for ${r}kb..."

    cool_dir="${imputed_matrix}/${r}kb_resolution/cool"
    if [[ ! -d "${cool_dir}" ]]; then
        echo "Error: Cool directory not found: ${cool_dir}" >&2
        exit 1
    fi

    find "${cool_dir}" -maxdepth 1 -type f -name "*.cool" -print \
        | LC_ALL=C sort > "${imputed_matrix}/${r}kb_resolution/filelist/valid_paths.txt"

    if [[ ! -s "${imputed_matrix}/${r}kb_resolution/filelist/valid_paths.txt" ]]; then
        echo "Error: No .cool files found in ${imputed_matrix}/${r}kb_resolution/cool/"
        exit 1
    fi

    awk -F'/' '{print $NF}' "${imputed_matrix}/${r}kb_resolution/filelist/valid_paths.txt" | sed 's/\.cool$//' > "${imputed_matrix}/${r}kb_resolution/filelist/ids.txt"

    paste "${imputed_matrix}/${r}kb_resolution/filelist/ids.txt" "${imputed_matrix}/${r}kb_resolution/filelist/valid_paths.txt" | sort -k1,1 > "${cell_table}"

    echo "Cell table updated. Total cells: $(wc -l < "${cell_table}")"

    mkdir -p "${imputed_matrix}/${r}kb_resolution/genescore"

    echo "Calculating gene score at ${r}kb..."
    hicluster gene-score \
        --cell_table_path "${cell_table}" \
        --gene_meta_path "${gene_meta}" \
        --resolution "${r}000" \
        --output_hdf_path "${imputed_matrix}/${r}kb_resolution/genescore/geneimputescore.hdf" \
        --chrom_size_path "${chrom_size}" \
        --cpu "${threads}" \
        --mode impute

    echo "Making AnnData at ${r}kb..."
    python "${SCRIPT_DIR}/01.pre-process/scripts/phc.make_adata.py" \
        -r "${r}" \
        --genome "${genome}" \
        "${chrom_size}" \
        "${imputed_matrix}/${r}kb_resolution/genescore/geneimputescore.hdf" \
        "${sample}"

    echo "Done with ${r}kb"
done
