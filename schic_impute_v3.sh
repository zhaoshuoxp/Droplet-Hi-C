#!/usr/bin/env bash

set -Eeuo pipefail

#######
current=$(pwd)
aln_dir=${current}/mapping
map_dir=${current}/splitPairs
in_dir=${current}/hicluster
raw_matrix=${in_dir}/raw_matrix
cell_matrix=${in_dir}/cell_matrix
imputed_matrix=${in_dir}/imputed_matrix

genome="mm10"
threads=24
res_array=(10)
s=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -s|--sample) s="$2"; shift ;;
        -t|--threads) threads="$2"; shift ;;
        -g|--genome) genome="$2"; shift ;;
        -r|--res) IFS=',' read -ra res_array <<< "$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$s" ]]; then
    echo "Sample prefix (-s) is required!"
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

#######

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOME_DIR="${DROPLET_HIC_GENOME_DIR:-/home/quanyiz/genome}"

mkdir -p "${map_dir}" "${cell_matrix}" "${raw_matrix}" "${imputed_matrix}"

case "${genome}" in
    mm10|hg38)
        chrom_size="${GENOME_DIR}/${genome}/${genome}_noYM.chrom.sizes"
        bl="${GENOME_DIR}/${genome}/${genome}-blacklist.v2.bed"
        ;;
    *)
        echo "Unsupported genome: ${genome}"
        exit 1
        ;;
esac

for required_file in "${chrom_size}" "${bl}" "${s}_cutoff.cells.txt" \
    "${aln_dir}/${s}_${genome}.sc.pairs.gz"; do
    if [[ ! -r "${required_file}" ]]; then
        echo "Required input is not readable: ${required_file}" >&2
        exit 1
    fi
done

cd "${current}"

echo -e "barcode\tlibrary\tsample\tcluster" > cluster.meta
awk -v SAMPLE="${s}" -v GENOME="${genome}" \
    'NF {print $1"\t"SAMPLE"\t"SAMPLE"_"GENOME".sc.pairs.gz\t"$1}' \
    "${s}_cutoff.cells.txt" >> cluster.meta
ulimit -n 10000 2>/dev/null || true
python "${SCRIPT_DIR}/01.pre-process/scripts/phc.sc_splitPairs_v3.py" \
    --indir "${aln_dir}" \
    --cluster cluster.meta \
    --outdir "${map_dir}" \
    --threads "${threads}" \
    --max 2048

awk -F '\t' -v OFS='\t' -v OUTDIR="${map_dir}" \
    'NR > 1 && !seen[$4]++ {print $4, OUTDIR "/" $4 ".pairs"}' \
    cluster.meta | LC_ALL=C sort -k1,1 > contact_table.tsv

if [[ ! -s contact_table.tsv ]]; then
    echo "No cells were selected in ${s}_cutoff.cells.txt." >&2
    exit 1
fi

while IFS=$'\t' read -r cell_id pairs_path; do
    if [[ ! -r "${pairs_path}" ]]; then
        echo "Split pairs file is missing for cell ${cell_id}: ${pairs_path}" >&2
        exit 1
    fi
done < contact_table.tsv

cell_table=${current}/contact_table.tsv

process_cell() {
    local cell_id=$1
    local pairs_path=$2
    local r c pad out_dir output_file chrom_dir
    echo "[$(date)] Start cell: ${cell_id}"

    printf '%s\t%s\n' "${cell_id}" "${pairs_path}" > "${raw_matrix}/${cell_id}.txt"

    hicluster filter-contact --output_dir "${raw_matrix}" --blacklist_1d_path "${bl}" \
        --chr1 1 --pos1 2 --chr2 3 --pos2 4 \
        --contact_table "${raw_matrix}/${cell_id}.txt" \
        --chrom_size_path "${chrom_size}" \
        --not_remove_duplicates

    for r in "${res_array[@]}"; do
        mkdir -p "${cell_matrix}/${r}kb_resolution" "${imputed_matrix}/${r}kb_resolution/cool"

        hicluster generatematrix-cell --infile "${raw_matrix}/${cell_id}.contact.rmbkl.tsv.gz" \
            --outdir "${cell_matrix}/${r}kb_resolution/" \
            --chrom_file "${chrom_size}" --res "${r}000" --cell "${cell_id}" \
            --chr1 1 --pos1 2 --chr2 3 --pos2 4

        pad=2
        if [[ "${r}" == "100" ]]; then
            pad=1
        fi

        shopt -s nullglob
        for chrom_dir in "${cell_matrix}/${r}kb_resolution"/chr*; do
            [[ ! -d "${chrom_dir}" ]] && continue
            c="${chrom_dir##*/}"
            [[ "${c}" == "chrM" || "${c}" == "chrY" ]] && continue

            mkdir -p "${imputed_matrix}/${r}kb_resolution/${c}"

            hicluster impute-cell --indir "${chrom_dir}/" \
                --outdir "${imputed_matrix}/${r}kb_resolution/${c}/" \
                --cell "${cell_id}" --chrom "${c}" --res "${r}000" --chrom_file "${chrom_size}" \
                --output_dist 500000000 --window_size 500000000 --step_size 500000000 \
                --pad "${pad}" --std 1 --rp 0.5 --tol 0.01

            out_dir="${imputed_matrix}/${r}kb_resolution/${c}/"
            for output_file in "${out_dir}"*.hdf5; do
                ln -sf "${output_file}" "${output_file%?}"
            done
        done
        shopt -u nullglob

        hic-internal aggregate-chromosomes \
            --chrom_size_path "${chrom_size}" \
            --resolution "${r}000" \
            --input_dir "${imputed_matrix}/${r}kb_resolution/" \
            --output_path "${imputed_matrix}/${r}kb_resolution/cool/${cell_id}.cool" \
            --chrom_wildcard "{chrom}/${cell_id}_{chrom}_pad${pad}_std1.0_rp0.5_sqrtvc.hdf" \
            --csr
    done

    echo "[$(date)] Finish cell: ${cell_id}"
}

pids=()
run_pool() {
    while (( $(jobs -rp | wc -l) >= threads )); do
        sleep 1
    done
    process_cell "$@" &
    pids+=("$!")
}

while IFS=$'\t' read -r cell_id pairs_path; do
    run_pool "${cell_id}" "${pairs_path}"
done < "${cell_table}"

pool_failed=0
for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
        pool_failed=1
    fi
done
if (( pool_failed != 0 )); then
    echo "One or more cell imputation jobs failed." >&2
    exit 1
fi
echo "✅ All jobs finished."
