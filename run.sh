#!/bin/bash

usage() {
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
    exit 1
}

if [ $# -eq 0 ]; then
    usage
fi


current=$(pwd)
threads=24
res="10,50,100"
steps="all" 
mode="atac"

while getopts "s:g:m:t:r:p:h" opt; do
  case ${opt} in
    s) s=${OPTARG} ;;              # Sample name
    g) genome=${OPTARG} ;;         # Genome (e.g., mm10, hg38)
    m) mode=${OPTARG} ;;           # Mode (e.g., atac, arc) 
    t) threads=${OPTARG} ;;        # Threads
    r) res=${OPTARG} ;;            # Resolution list i.e. 10,50,100
    p) steps=${OPTARG} ;;          # Steps to run
    h) usage ;;                    # Show help
    *) usage ;;
  esac
done


requires=("hictools" "trim_galore" "bwa" "samtools" "pairtools" "bgzip" "pairix" "cooler" "hicluster" "hic-internal")
for i in ${requires[@]};do
	which $i &>/dev/null || { echo "Error: $i not found"; exit 1; }
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $current

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

# Step 1: preproc_paired_hic
if run_step 1; then
    echo ">>> Running Step 1: preproc_paired_hic..."
    $SCRIPT_DIR/preproc_paired_hic_v3.sh -s $s -g $genome -m $mode -t $threads
fi

# Step 2: cell_filtering
if run_step 2; then
    echo ">>> Running Step 2: cell_filtering..."
    $SCRIPT_DIR/01.pre-process/scripts/phc.cell_filtering.R $s
fi

# Step 3: schic_impute
if run_step 3; then
    echo ">>> Running Step 3: schic_impute..."
    $SCRIPT_DIR/schic_impute_v3.sh -s $s -t $threads -g $genome -r $res
fi
 
# Step 4: schic_gene_score
if run_step 4; then
    echo ">>> Running Step 4: schic_gene_score..."
    $SCRIPT_DIR/schic_gene_score_v3.sh -s $s -t $threads -g $genome -r $res
fi

# Step 5: schic_domain
if run_step 5; then
    echo ">>> Running Step 5: schic_domain..."
    $SCRIPT_DIR/schic_domain_v3.sh -s $s -t $threads -r $res
fi

# Step 6: schic_loops
if run_step 6; then
    echo ">>> Running Step 6: schic_loops..."
    $SCRIPT_DIR/schic_loops_v3.sh -t $threads -g $genome -r $res
fi

echo ">>> All specified tasks completed successfully!"