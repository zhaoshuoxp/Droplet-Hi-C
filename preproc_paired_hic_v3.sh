#!/usr/bin/env bash

set -Eeuo pipefail

current=$(pwd)
s=""
genome="mm10"
mode="atac"
threads=24

# Parse command-line arguments
while getopts "s:g:m:o:t:" opt; do
  case ${opt} in
    s) s=${OPTARG} ;;  # Sample name
    g) genome=${OPTARG} ;;  # Genome (e.g., mm10, hg38)
    m) mode=${OPTARG} ;;  # Mode (e.g., atac, arc)
    o) current=${OPTARG} ;;  # Output directory
    t) threads=${OPTARG} ;;
    *) echo "Usage: $0 [-s sample] [-g genome] [-m mode] [-o output] [-t 24]"; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOME_DIR="${DROPLET_HIC_GENOME_DIR:-/home/quanyiz/genome}"

################################
fastq_dir="${current}/fastq/"
trim_dir="${current}/trimmed/"
map_dir="${current}/mapping/"
mtx_dir="${current}/matrices/"
script_dir="${SCRIPT_DIR}/01.pre-process/scripts"

mm10_bwa="${GENOME_DIR}/mm10/BWAindex/mm10bwa"
hg38_bwa="${GENOME_DIR}/hg38/BWAindex/hg38bwa"
mm10_chrsize="${GENOME_DIR}/mm10/mm10_noYM.chrom.sizes"
hg38_chrsize="${GENOME_DIR}/hg38/hg38_noYM.chrom.sizes"
atac_10X="${SCRIPT_DIR}/01.pre-process/10xBC_index/737K-cratac-v1"
arc_10X="${SCRIPT_DIR}/01.pre-process/10xBC_index/737K-arc-v1"
##################################

if [[ -z "${s}" ]]; then
    echo "Error: sample name (-s) is required." >&2
    exit 1
fi
if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: threads must be a positive integer." >&2
    exit 1
fi

case "${genome}" in
    mm10) ref="${mm10_bwa}"; chrsize="${mm10_chrsize}" ;;
    hg38) ref="${hg38_bwa}"; chrsize="${hg38_chrsize}" ;;
    *) echo "Error: unsupported genome '${genome}' (expected mm10 or hg38)." >&2; exit 1 ;;
esac

case "${mode}" in
    atac) ref_10X="${atac_10X}" ;;
    arc) ref_10X="${arc_10X}" ;;
    *) echo "Error: unsupported mode '${mode}' (expected atac or arc)." >&2; exit 1 ;;
esac

hictools_bin="${DROPLET_HIC_HICTOOLS:-${SCRIPT_DIR}/hictools_2/hictools}"
if [[ ! -x "${hictools_bin}" ]]; then
    if ! hictools_bin="$(command -v hictools)"; then
        echo "Error: hictools executable not found." >&2
        exit 1
    fi
fi

for required_file in \
    "${ref}.bwt" \
    "${chrsize}" \
    "${ref_10X}.fa"; do
    if [[ ! -r "${required_file}" ]]; then
        echo "Error: required input is not readable: ${required_file}" >&2
        exit 1
    fi
done

for read_number in 1 2 3; do
    if [[ ! -r "${fastq_dir}/${s}_R${read_number}.fq.gz" \
          && ! -r "${fastq_dir}/${s}_R${read_number}.fastq.gz" ]]; then
        echo "Error: cannot find readable R${read_number} FASTQ for sample '${s}' in ${fastq_dir}" >&2
        exit 1
    fi
done

mkdir -p "${trim_dir}" "${map_dir}" "${mtx_dir}"

echo "process scHiC fastq... mode: "${mode}
### NEW hictools end-to-end 1 step processing
"${hictools_bin}" end_to_end "${mode}" "${fastq_dir}/${s}" "${ref_10X}.fa" "${threads}"

trim_galore -q 20 -j "${threads}" --paired \
    "${fastq_dir}/${s}_R1_BC_cov.fq.gz" \
    "${fastq_dir}/${s}_R3_BC_cov.fq.gz" \
    -o "${trim_dir}"

# ### mapping
(bwa mem -SP5M -T0 -t "${threads}" "${ref}" \
    "${trim_dir}/${s}_R1_BC_cov_val_1.fq.gz" \
    "${trim_dir}/${s}_R3_BC_cov_val_2.fq.gz" \
    | samtools view -bhS - > "${map_dir}/${s}_${genome}.bam") \
    2>"${map_dir}/${s}_${genome}.log"

### identify valid pair interactions using pairtools
### all possible complex ligations are preserved (--walks-policy all)
samtools view -h -@ "${threads}" "${map_dir}/${s}_${genome}.bam" | \
pairtools parse \
    --min-mapq 40 \
    --walks-policy all \
    --nproc-in "${threads}" \
    --nproc-out "${threads}" \
    --max-inter-align-gap 30 \
    --chroms-path "${chrsize}" \
    --assembly "${genome}" \
    --output-stats "${map_dir}/${s}_${genome}.pairparse.txt" | \
pairtools sort --nproc "${threads}" --tmpdir="${map_dir}" | \
awk -v FS="\t" -v OFS="\t" '
    BEGIN { corrupt = 0 }
    /^#columns/ { print $0 " barcode1 barcode2"; next }
    /^#/ { print; next }
    {
        if ($8 == "XX" || NF != 10) {
            corrupt++;
            next;
        }
        split($9, a, /[:^]+/); bc1 = a[1];
        split($10, b, /[:^]+/); bc2 = b[1];

        if (length(bc1) != 16 || length(bc2) != 16) {
            print "Error: Barcodes at head of readname is not 16bp at line " NR > "/dev/stderr";
            exit 1;
        }
        print $0, bc1, bc2
    }
    END {
        print "Pairs removed due to unmapped (type: XX): " corrupt > "/dev/stderr"
    }
' | \
pairtools dedup \
    --nproc-in "${threads}" \
    --nproc-out "${threads}" \
    --extra-col-pair "barcode1" "barcode2" \
    --mark-dups \
    --output-stats "${map_dir}/${s}_${genome}.sc.pairdedup.txt" | \
pairtools split \
    --nproc-in "${threads}" \
    --nproc-out "${threads}" \
    --output-pairs "${map_dir}/${s}_${genome}.sc.pairs" \
    --output-sam - | \
samtools view -b -@ "${threads}" | \
samtools sort -T "${map_dir}/samtools-sort" -@ "${threads}" -m 4G \
    -o "${map_dir}/${s}_${genome}.sc.pairtools.bam"

Rscript "${script_dir}/phc.summarize_pairs_lec.R" "${map_dir}/${s}_${genome}.sc.pairdedup.txt"
Rscript "${script_dir}/phc.plot_fragment.R" "${map_dir}/${s}_${genome}.sc.pairs"

bgzip -f "${map_dir}/${s}_${genome}.sc.pairs"
pairix -f "${map_dir}/${s}_${genome}.sc.pairs.gz"

### count valid pair per barcode (valid reads type?)
python "${script_dir}/phc.count_pairs_sc.py" \
    --input "${map_dir}/${s}_${genome}.sc.pairs.gz" \
    --output "${map_dir}/${s}_${genome}.PairCount"

### binning
bs=5000
cooler cload pairix "${chrsize}:${bs}" \
    "${map_dir}/${s}_${genome}.sc.pairs.gz" \
    "${mtx_dir}/${s}_${genome}_${bs}.cool"

### Balancing
cooler zoomify --balance -p "${threads}" \
    -o "${mtx_dir}/${s}_${genome}.mcool" \
    -r 5000,10000,25000,50000,100000,250000,500000,1000000,2500000 \
    "${mtx_dir}/${s}_${genome}_${bs}.cool"
