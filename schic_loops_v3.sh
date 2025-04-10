#!/bin/bash

genome="mm10"
threads=24
res_array=(10)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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
    chrom_size=$SCRIPT_DIR"/01.pre-process/supp/mm10.chrom.sizes"
elif [[ "${genome}" == "hg38" ]]
then
    chrom_size=$SCRIPT_DIR"/01.pre-process/supp/hg38.chrom.sizes"
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
            [[ "$c" == "chr*" ]] && break
            [[ "$c" == "chrM" || "$c" == "chrY" ]] && continue
            echo "Processing resolution ${r}kb for cell ${s} (chromosome ${c})"
            hicluster loop-bkg-cell --indir "${imputed_matrix}/${r}kb_resolution/" \
            --cell ${s} \
            --chrom ${c} \
            --impute_mode pad2_std1.0_rp0.5_sqrtvc \
            --res "${r}000" \
            --pad 2 \
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
        ls "${imputed_matrix}/${r}kb_resolution/${c}" | grep "pad2_std1.0_rp0.5_sqrtvc.hdf5" - | grep -f ${meta} - | sed 's/.hdf5//g' - > "${imputed_matrix}/${r}kb_resolution/filelist/imputematrices_${c}_1.txt"
        awk -v path="${imputed_matrix}/${r}kb_resolution/${c}/" '{print path $0}' "${imputed_matrix}/${r}kb_resolution/filelist/imputematrices_${c}_1.txt" > "${imputed_matrix}/${r}kb_resolution/filelist/imputematrices_${c}.txt"
        
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
    hicluster loop-mergechr --inprefix "${imputed_matrix}/${r}kb_resolution/matrices/imputematrices" \
        --outprefix "${imputed_matrix}/${r}kb_resolution/loop/imputematrices" \
        --chrom_file "${chrom_size}"
done

echo "✅ Loop merge finished for all resolutions."
