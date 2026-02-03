#!/bin/bash

genome="mm10"
threads=24
res_array=(50)

while getopts "g:r:t:" opt; do
  case ${opt} in
    g) genome=${OPTARG} ;;
    r) IFS=',' read -ra res_array <<< "${OPTARG}" ;;  
    t) threads=${OPTARG} ;; 
    *) echo "Usage: $0 [-g genome] [-r res1,res2,...] [-t threads]"; exit 1 ;;
  esac
done

current=$(pwd)
imputed_matrix="${current}/hicluster/imputed_matrix"
cut -f1 ${current}/contact_table.tsv > cell.id
meta="${current}/cell.id"  

if [[ "${genome}" == "mm10" ]]
then
    chrom_size="/home/quanyiz/genome/mm10/mm10_noYM.chrom.sizes"
elif [[ "${genome}" == "hg38" ]]
then
    chrom_size="/home/quanyiz/genome/hg38/hg38_noYM.chrom.sizes"
fi

for r in "${res_array[@]}"; do
    mkdir -p "${imputed_matrix}/${r}kb_resolution/matrices/" "${imputed_matrix}/${r}kb_resolution/loop/"
done

chroms=$(cut -f1 $chrom_size)

start=1
end=$(wc -l $meta | awk '{print $1}')

loop_bkg_cell() {
    i=$1
    s=$(awk -v iline=${i} 'NR==iline {print $1}' ${meta})
    echo "[`date`] Start cell: $s (line $i)"
    for r in "${res_array[@]}"; do
        for c in $chroms; do
            [[ "$c" != "chr"* ]] && continue
            [[ "$c" == "chrM" || "$c" == "chrY" ]] && continue

            file_pattern="${imputed_matrix}/${r}kb_resolution/${c}/${s}_${c}_pad?_std1.0_rp0.5_sqrtvc.hdf5"

            files=( $file_pattern )
            input_file="${files[0]}"

            if [[ ! -f "$input_file" ]]; then
                continue
            fi

            if [[ "$input_file" == *"_pad1_"* ]]; then
                pad_val=1
            else
                pad_val=2
            fi

            echo "Processing resolution ${r}kb for cell ${s} (chromosome ${c}) using pad${pad_val}"
            hicluster loop-bkg-cell --indir "${imputed_matrix}/${r}kb_resolution/" \
            --cell ${s} \
            --chrom ${c} \
            --impute_mode "pad${pad_val}_std1.0_rp0.5_sqrtvc" \
            --res "${r}000" \
            --pad ${pad_val} \
            --dist 10000000
        done
    done
    echo "[`date`] Finish cell: $s"
}

run_pool_bkg() {
    while (( $(jobs -rp | wc -l) >= threads )); do
        sleep 1
    done
    loop_bkg_cell "$@" &
}

for ((i=start; i<=end; i++)); do
    run_pool_bkg "$i"
done

wait
echo "✅ Background processing finished."

loop_sumcell_cell() {
    c=$1
    for r in "${res_array[@]}"; do

        echo "Generating file list for ${c}..."
        find "${imputed_matrix}/${r}kb_resolution/${c}" -name "*_dist_trim.E.npz" \
        | sed "s/_dist_trim.E.npz//g" \
        | grep -F -f "${meta}" \
        > "${imputed_matrix}/${r}kb_resolution/filelist/imputematrices_${c}.txt"

        if [[ ! -s "${imputed_matrix}/${r}kb_resolution/filelist/imputematrices_${c}.txt" ]]; then
            echo "Warning: No valid matrix files found for ${c}. Skipping."
            continue
        fi
        
        echo "Summing loop for resolution ${r}kb on chromosome ${c}..."
        hicluster loop-sumcell-chr --cell_list "${imputed_matrix}/${r}kb_resolution/filelist/imputematrices_${c}.txt" \
            --outprefix "${imputed_matrix}/${r}kb_resolution/matrices/imputematrices_${c}" \
            --res "${r}000"
    done
}

run_pool_sumcell() {
    while (( $(jobs -rp | wc -l) >= threads )); do
        sleep 1
    done
    loop_sumcell_cell "$@" &
}

for c in $chroms; do
    run_pool_sumcell "$c"
done

wait
echo "✅ Cell summing finished."

for r in "${res_array[@]}"; do
    echo "Merging results for resolution ${r}kb..."
    dist_val=$((r * 2000))
    
    hicluster loop-mergechr \
        --inprefix "${imputed_matrix}/${r}kb_resolution/matrices/imputematrices" \
        --outprefix "${imputed_matrix}/${r}kb_resolution/loop/imputematrices" \
        --chrom_file "${chrom_size}" \
        --res "${r}000" \
        --dist_thres "${dist_val}" \
        --fdr_thres 0.1
done

echo "✅ Loop merge finished for all resolutions."
