#!/usr/bin/env bash

set -Eeuo pipefail

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
if [[ ! -r "${current}/contact_table.tsv" ]]; then
    echo "Contact table not found: ${current}/contact_table.tsv" >&2
    exit 1
fi
cut -f1 "${current}/contact_table.tsv" > cell.id
meta="${current}/cell.id"

GENOME_DIR="${DROPLET_HIC_GENOME_DIR:-/home/quanyiz/genome}"
if [[ "${genome}" != "mm10" && "${genome}" != "hg38" ]]; then
    echo "Unsupported genome: ${genome}" >&2
    exit 1
fi
chrom_size="${GENOME_DIR}/${genome}/${genome}_noYM.chrom.sizes"
if [[ ! -r "${chrom_size}" ]]; then
    echo "Chromosome sizes file is not readable: ${chrom_size}" >&2
    exit 1
fi

mapfile -t cells < "${meta}"
if (( ${#cells[@]} == 0 )); then
    echo "No cells found in ${current}/contact_table.tsv" >&2
    exit 1
fi

for r in "${res_array[@]}"; do
    mkdir -p \
        "${imputed_matrix}/${r}kb_resolution/filelist/" \
        "${imputed_matrix}/${r}kb_resolution/matrices/" \
        "${imputed_matrix}/${r}kb_resolution/loop/"
done

mapfile -t chroms < <(cut -f1 "${chrom_size}")

loop_bkg_cell() {
    local s=$1
    local r c input_file pad_val
    echo "[$(date)] Start cell: ${s}"
    for r in "${res_array[@]}"; do
        pad_val=2
        if [[ "${r}" == "100" ]]; then
            pad_val=1
        fi
        for c in "${chroms[@]}"; do
            [[ "${c}" != chr* ]] && continue
            [[ "${c}" == "chrM" || "${c}" == "chrY" ]] && continue

            input_file="${imputed_matrix}/${r}kb_resolution/${c}/${s}_${c}_pad${pad_val}_std1.0_rp0.5_sqrtvc.hdf5"

            if [[ ! -f "${input_file}" ]]; then
                continue
            fi

            echo "Processing resolution ${r}kb for cell ${s} (chromosome ${c}) using pad${pad_val}"
            hicluster loop-bkg-cell --indir "${imputed_matrix}/${r}kb_resolution/" \
                --cell "${s}" \
                --chrom "${c}" \
                --impute_mode "pad${pad_val}_std1.0_rp0.5_sqrtvc" \
                --res "${r}000" \
                --pad "${pad_val}" \
                --dist 10000000
        done
    done
    echo "[$(date)] Finish cell: ${s}"
}

pids=()
run_pool_bkg() {
    while (( $(jobs -rp | wc -l) >= threads )); do
        sleep 1
    done
    loop_bkg_cell "$@" &
    pids+=("$!")
}

for cell_id in "${cells[@]}"; do
    run_pool_bkg "${cell_id}"
done

pool_failed=0
for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
        pool_failed=1
    fi
done
if (( pool_failed != 0 )); then
    echo "One or more loop background jobs failed." >&2
    exit 1
fi
echo "✅ Background processing finished."

loop_sumcell_cell() {
    local c=$1
    local r chrom_dir output_list
    for r in "${res_array[@]}"; do

        echo "Generating file list for ${c}..."
        chrom_dir="${imputed_matrix}/${r}kb_resolution/${c}"
        output_list="${imputed_matrix}/${r}kb_resolution/filelist/imputematrices_${c}.txt"
        if [[ ! -d "${chrom_dir}" ]]; then
            : > "${output_list}"
        else
            find "${chrom_dir}" -maxdepth 1 -type f -name "*_dist_trim.E.npz" -print \
                | sed 's/_dist_trim\.E\.npz$//' \
                | awk -v CHROM="${c}" '
                    NR == FNR { keep[$1] = 1; next }
                    {
                        path = $0
                        name = path
                        sub(/^.*\//, "", name)
                        sub("_" CHROM "$", "", name)
                        if (name in keep) print path
                    }
                ' "${meta}" - \
                | LC_ALL=C sort > "${output_list}"
        fi

        if [[ ! -s "${output_list}" ]]; then
            echo "Warning: No valid matrix files found for ${c}. Skipping."
            continue
        fi

        echo "Summing loop for resolution ${r}kb on chromosome ${c}..."
        hicluster loop-sumcell-chr --cell_list "${output_list}" \
            --outprefix "${imputed_matrix}/${r}kb_resolution/matrices/imputematrices_${c}" \
            --res "${r}000"
    done
}

pids=()
run_pool_sumcell() {
    while (( $(jobs -rp | wc -l) >= threads )); do
        sleep 1
    done
    loop_sumcell_cell "$@" &
    pids+=("$!")
}

for c in "${chroms[@]}"; do
    [[ "${c}" != chr* ]] && continue
    [[ "${c}" == "chrM" || "${c}" == "chrY" ]] && continue
    run_pool_sumcell "${c}"
done

pool_failed=0
for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
        pool_failed=1
    fi
done
if (( pool_failed != 0 )); then
    echo "One or more chromosome summing jobs failed." >&2
    exit 1
fi
echo "✅ Cell summing finished."

for r in "${res_array[@]}"; do
    echo "Merging results for resolution ${r}kb..."
    dist_val=$((10#${r} * 2000))

    hicluster loop-mergechr \
        --inprefix "${imputed_matrix}/${r}kb_resolution/matrices/imputematrices" \
        --outprefix "${imputed_matrix}/${r}kb_resolution/loop/imputematrices" \
        --chrom_file "${chrom_size}" \
        --res "${r}000" \
        --dist_thres "${dist_val}" \
        --fdr_thres 0.1
done

echo "✅ Loop merge finished for all resolutions."
