#!/bin/bash
# check dependences
requires=("hictools" "trim_galore" "bwa" "samtools" "pairtools" "bgzip" "pairix" "cooler" "hicluster" "hic-internal")
for i in ${requires[@]};do
	which $i &>/dev/null || { echo $i not found; exit 1; }
done

current=$(pwd)
threads=24
res=10
# Parse command-line arguments
while getopts "s:g:m:t:r:" opt; do
  case ${opt} in
    s) s=${OPTARG} ;;              # Sample name
    g) genome=${OPTARG} ;;         # Genome (e.g., mm10, hg38)
    m) mode=${OPTARG} ;;           # Mode (e.g., atac, arc) 
    t) threads=${OPTARG} ;;        # Threads
    r) res=${OPTARG} ;;            # Resolution list i.e. 10,50,100
    *) echo "Usage: $0 [-s sample] [-g genome] [-m mode] [-t threads] [-r res1,res2,...]"; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd $current

$SCRIPT_DIR/preproc_paired_hic_v3.sh -s $s -g $genome -m $mode -t $threads

$SCRIPT_DIR/01.pre-process/scripts/phc.cell_filtering.R $s

$SCRIPT_DIR/schic_impute_v3.sh -s $s -t $threads -g $genome -r $res
 
$SCRIPT_DIR/schic_gene_score_v3.sh -s $s -t $threads -g $genome -r $res

$SCRIPT_DIR/schic_domain_v3.sh -s $s -t $threads -r $res

$SCRIPT_DIR/schic_loops_v3.sh -t $threads -g $genome -r $res
