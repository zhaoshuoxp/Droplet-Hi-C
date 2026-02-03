#!/bin/bash

genome="mm10"
threads=24
res_array=(10)

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

if [[ "${genome}" == "mm10" ]]; then
    gene_meta="/nfs/baldar/quanyiz/app/Droplet-Hi-C/01.pre-process/bed/mm10_gene_anno.bed"
    chrom_size="/home/quanyiz/genome/mm10/mm10_noYM.chrom.sizes"
elif [[ "${genome}" == "hg38" ]]; then
    gene_meta="/nfs/baldar/quanyiz/app/Droplet-Hi-C/01.pre-process/bed/hg38_gene_anno.bed"
    chrom_size="/home/quanyiz/genome/hg38/hg38_noYM.chrom.sizes"
else
    echo "Unsupported genome: ${genome}"
    exit 1
fi

current=$(pwd)
imputed_matrix="${current}/hicluster/imputed_matrix"

for r in "${res_array[@]}"; do
    echo "Processing resolution ${r}kb..."

    mkdir -p "${imputed_matrix}/${r}kb_resolution/filelist"
    cell_table="${imputed_matrix}/${r}kb_resolution/filelist/cell_table.tsv"

    echo "Generating/Updating cell table for ${r}kb..."

    find "${imputed_matrix}/${r}kb_resolution/cool" -name "*.cool" > "${imputed_matrix}/${r}kb_resolution/filelist/valid_paths.txt"

    if [[ ! -s "${imputed_matrix}/${r}kb_resolution/filelist/valid_paths.txt" ]]; then
        echo "Error: No .cool files found in ${imputed_matrix}/${r}kb_resolution/cool/"
        continue
    fi

    awk -F'/' '{print $NF}' "${imputed_matrix}/${r}kb_resolution/filelist/valid_paths.txt" | sed 's/\.cool$//' > "${imputed_matrix}/${r}kb_resolution/filelist/ids.txt"

    paste "${imputed_matrix}/${r}kb_resolution/filelist/ids.txt" "${imputed_matrix}/${r}kb_resolution/filelist/valid_paths.txt" | sort -k1,1 > "${cell_table}"

    echo "Cell table updated. Total cells: $(wc -l < ${cell_table})"

    mkdir -p "${imputed_matrix}/${r}kb_resolution/genescore"

    echo "Calculating gene score at ${r}kb..."
    hicluster gene-score \
        --cell_table_path "${cell_table}" \
        --gene_meta_path "${gene_meta}" \
        --resolution "${r}000" \
        --output_hdf_path "${imputed_matrix}/${r}kb_resolution/genescore/geneimputescore.hdf" \
        --chrom_size_path "${chrom_size}" \
        --cpu "$threads" \
        --mode impute
        
    echo "Making AnnData at ${r}kb..."
    /nfs/baldar/quanyiz/app/Droplet-Hi-C/01.pre-process/scripts/phc.make_adata.py \
    -r $r \
    --genome $genome \
    $chrom_size \
    ${imputed_matrix}/${r}kb_resolution/genescore/geneimputescore.hdf \
    $sample
    
    echo "Done with ${r}kb"
done