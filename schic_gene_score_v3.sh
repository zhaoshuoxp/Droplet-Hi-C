#!/bin/bash

genome="mm10"
threads=24
res_array=(10)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while getopts "s:g:r:t:" opt; do
  case ${opt} in
    s) sample=${OPTARG} ;;
    g) genome=${OPTARG} ;;
    r) IFS=',' read -ra res_array <<< "${OPTARG}" ;;
    t) threads=${OPTARG} ;;
    *) echo "Usage: $0 [-g genome] [-r res1,res2,...] [-t threads]"; exit 1 ;;
  esac
done

if [ -z "${sample}" ]; then
  echo "Sample name (-s) is required."
  exit 1
fi

if [[ "${genome}" == "mm10" ]]; then
    gene_meta=$SCRIPT_DIR"/01.pre-process/bed/mm10_gene_anno.bed"
    chrom_size=$SCRIPT_DIR"/01.pre-process/supp/mm10.chrom.sizes"
elif [[ "${genome}" == "hg38" ]]; then
    gene_meta=$SCRIPT_DIR"/01.pre-process/bed/hg38_gene_anno.bed"
    chrom_size=$SCRIPT_DIR"/01.pre-process/supp/hg38.chrom.sizes"
else
    echo "Unsupported genome: ${genome}"
    exit 1
fi

current=$(pwd)
imputed_matrix="${current}/hicluster/imputed_matrix"

for r in "${res_array[@]}"; do
    echo "Processing resolution ${r}kb..."

    cell_table="${imputed_matrix}/${r}kb_resolution/filelist/cell_table.tsv"

    if [[ -f ${cell_table} ]]; then
        echo "Cell table exists for ${r}kb"
    else
        mkdir -p "${imputed_matrix}/${r}kb_resolution/filelist"
        ls "${imputed_matrix}/${r}kb_resolution/cool/" | grep '.cool' | \
        awk -v p="${imputed_matrix}/${r}kb_resolution/cool" '{printf("%s/%s\n", p, $0)}' > "${imputed_matrix}/${r}kb_resolution/filelist/cell.txt"

        paste <(awk -F'/' '{print $NF}' "${imputed_matrix}/${r}kb_resolution/filelist/cell.txt" | cut -d. -f1) \
              "${imputed_matrix}/${r}kb_resolution/filelist/cell.txt" | sort -k1,1 \
              > "${imputed_matrix}/${r}kb_resolution/filelist/cell.txt.tmp"
        
        mv "${imputed_matrix}/${r}kb_resolution/filelist/cell.txt.tmp" "${cell_table}"
    fi

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
    $SCRIPT_DIR/01.pre-process/scripts/phc.make_adata.py \
    -r $r \
    --genome $genome \
    $chrom_size \
    ${imputed_matrix}/${r}kb_resolution/genescore/geneimputescore.hdf \
    $sample
    
    echo "Done with ${r}kb"
done
