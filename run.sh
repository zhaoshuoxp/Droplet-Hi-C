#!/usr/bin/env bash

set -Eeuo pipefail

usage() {
    local exit_code="${1:-1}"
    echo "Usage: $0 [-s sample] [-g genome] [-m mode] [-t threads] [-r res1,res2,...] [-p steps] [-h]"
    echo ""
    echo "Options:"
    echo "  -s  Sample name"
    echo "  -g  Genome (e.g., mm10, hg38)"
    echo "  -m  Mode (e.g., atac, arc, default: atac)"
    echo "  -t  Threads (default: 24)"
    echo "  -r  Resolution list i.e. 10 (default: 10,50,100)"
    echo "  -p  Steps to run, comma-separated (default: all)"
    echo "        1 : preproc_paired_hic"
    echo "        2 : cell_filtering"
    echo "        3 : schic_impute"
    echo "        4 : schic_gene_score"
    echo "        5 : schic_domain"
    echo "        6 : schic_loops"
    echo "        Example: -p 1,2,3 or -p 4,6"
    echo "  -h  Show this help message and exit"
    exit "${exit_code}"
}

if [ $# -eq 0 ]; then
    usage
fi


current=$(pwd)
threads=24
res="10,50,100"
steps="all"
mode="atac"
s=""
genome=""

while getopts "s:g:m:t:r:p:h" opt; do
  case ${opt} in
    s) s=${OPTARG} ;;              # Sample name
    g) genome=${OPTARG} ;;         # Genome (e.g., mm10, hg38)
    m) mode=${OPTARG} ;;           # Mode (e.g., atac, arc)
    t) threads=${OPTARG} ;;        # Threads
    r) res=${OPTARG} ;;            # Resolution list i.e. 10,50,100
    p) steps=${OPTARG} ;;          # Steps to run
    h) usage 0 ;;                  # Show help
    *) usage ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${current}"

run_step() {
    local step_num=$1
    if [[ "$steps" == "all" ]]; then
        return 0
    fi

    if [[ ",${steps}," == *",${step_num},"* ]]; then
        return 0
    else
        return 1
    fi
}

sample_required=0
for step_num in 1 2 3 4 5; do
    if run_step "${step_num}"; then
        sample_required=1
        break
    fi
done

genome_required=0
for step_num in 1 2 3 4 6; do
    if run_step "${step_num}"; then
        genome_required=1
        break
    fi
done

if (( sample_required == 1 )) && [[ -z "${s}" ]]; then
    echo "Error: -s sample is required for the selected steps." >&2
    usage
fi

if (( genome_required == 1 )) && [[ -z "${genome}" ]]; then
    echo "Error: -g genome is required for the selected steps." >&2
    usage
fi

if (( genome_required == 1 )) && [[ "${genome}" != "mm10" && "${genome}" != "hg38" ]]; then
    echo "Error: unsupported genome '${genome}' (expected mm10 or hg38)." >&2
    exit 1
fi

if run_step 1 && [[ "${mode}" != "atac" && "${mode}" != "arc" ]]; then
    echo "Error: unsupported mode '${mode}' (expected atac or arc)." >&2
    exit 1
fi

if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: threads must be a positive integer." >&2
    exit 1
fi

if [[ ! "${res}" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ]]; then
    echo "Error: resolutions must be comma-separated positive integers." >&2
    exit 1
fi

if [[ "${steps}" != "all" && ! "${steps}" =~ ^[1-6](,[1-6])*$ ]]; then
    echo "Error: steps must be 'all' or a comma-separated list containing 1-6." >&2
    exit 1
fi

declare -A required_commands=()
add_requirement() {
    required_commands["$1"]=1
}

if run_step 1; then
    for command_name in trim_galore bwa samtools pairtools bgzip pairix cooler Rscript python; do
        add_requirement "${command_name}"
    done
    if [[ ! -x "${SCRIPT_DIR}/hictools_2/hictools" ]]; then
        add_requirement hictools
    fi
fi
if run_step 2; then
    add_requirement Rscript
fi
if run_step 3; then
    add_requirement python
    add_requirement hicluster
    add_requirement hic-internal
fi
if run_step 4; then
    add_requirement python
    add_requirement hicluster
fi
if run_step 5 || run_step 6; then
    add_requirement hicluster
fi

for command_name in "${!required_commands[@]}"; do
    command -v "${command_name}" >/dev/null 2>&1 || {
        echo "Error: ${command_name} not found" >&2
        exit 1
    }
done

# Step 1: preproc_paired_hic
if run_step 1; then
    echo ">>> Running Step 1: preproc_paired_hic..."
    "${SCRIPT_DIR}/preproc_paired_hic_v3.sh" -s "${s}" -g "${genome}" -m "${mode}" -t "${threads}"
fi

# Step 2: cell_filtering
if run_step 2; then
    echo ">>> Running Step 2: cell_filtering..."
    Rscript "${SCRIPT_DIR}/01.pre-process/scripts/phc.cell_filtering.R" "${s}" "${genome}"
fi

# Step 3: schic_impute
if run_step 3; then
    echo ">>> Running Step 3: schic_impute..."
    "${SCRIPT_DIR}/schic_impute_v3.sh" -s "${s}" -t "${threads}" -g "${genome}" -r "${res}"
fi

# Step 4: schic_gene_score
if run_step 4; then
    echo ">>> Running Step 4: schic_gene_score..."
    "${SCRIPT_DIR}/schic_gene_score_v3.sh" -s "${s}" -t "${threads}" -g "${genome}" -r "${res}"
fi

# Step 5: schic_domain
if run_step 5; then
    echo ">>> Running Step 5: schic_domain..."
    "${SCRIPT_DIR}/schic_domain_v3.sh" -s "${s}" -t "${threads}" -r "${res}"
fi

# Step 6: schic_loops
if run_step 6; then
    echo ">>> Running Step 6: schic_loops..."
    "${SCRIPT_DIR}/schic_loops_v3.sh" -t "${threads}" -g "${genome}" -r "${res}"
fi

echo ">>> All specified tasks completed successfully!"
