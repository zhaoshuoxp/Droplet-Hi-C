#!/bin/bash

current=$(pwd)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

################################
fastq_dir="${current}/fastq/"
trim_dir="${current}/trimmed/"
map_dir="${current}/mapping/"
old_map_dir="${current}/prev_mapping/"
mtx_dir="${current}/matrices/"
script_dir=$SCRIPT_DIR"/01.pre-process/scripts"

mm10_bwa="/home/quanyiz/genome/mm10/BWAindex/mm10bwa"
hg38_bwa="/home/quanyiz/genome/hg38/BWAindex/hg38bwa"
mm10_chrsize=$SCRIPT_DIR"/01.pre-process/supp/mm10.chrom.sizes"
hg38_chrsize=$SCRIPT_DIR"/01.pre-process/supp/hg38.chrom.sizes"
mm10_bl=$SCRIPT_DIR"/01.pre-process/supp/mm10-blacklist.v2.bed"
hg38_bl=$SCRIPT_DIR"/01.pre-process/supp/hg38-blacklist.v2.bed" 
atac_10X=$SCRIPT_DIR"/01.pre-process/10xBC_index/737K-cratac-v1"
arc_10X=$SCRIPT_DIR"/01.pre-process/10xBC_index/737K-arc-v1"
##################################

if [ -z ${genome+x} ]; then echo "genome is not specified. Default to mm10"; genome="mm10"; fi
if [ -z ${mode+x} ]; then echo "10X kit used is not specified. Default to atac"; mode="atac"; fi
if [ -z ${threads+x} ]; then echo "threads used is not specified. Default to 24"; threads=24; fi

if [ $genome == "mm10" ]; then ref=${mm10_bwa}; chrsize=${mm10_chrsize}; fi
if [ $genome == "hg38" ]; then ref=${hg38_bwa}; chrsize=${hg38_chrsize}; fi
if [ $mode == "arc" ]; then ref_10X=${arc_10X}; else ref_10X=${atac_10X}; fi

mkdir -p ${trim_dir}
mkdir -p ${map_dir}
mkdir -p ${mtx_dir}

echo "process scHiC fastq... mode: "${mode}
hictools combine_hic ${mode} ${fastq_dir}/${s}
# there is no --reorder in bowtie1, use single-thread, NO -p param!
(zcat ${fastq_dir}/${s}_R1_combined.fq.gz | bowtie ${ref_10X} - --nofw -m 1 -v 1 -S ${fastq_dir}/${s}_R1_BC.sam) 2>${fastq_dir}/${s}_R1.log
(zcat ${fastq_dir}/${s}_R3_combined.fq.gz | bowtie ${ref_10X} - --nofw -m 1 -v 1 -S ${fastq_dir}/${s}_R3_BC.sam) 2>${fastq_dir}/${s}_R3.log

hictools convert_hic2 ${fastq_dir}/${s}_R1_BC.sam
hictools convert_hic2 ${fastq_dir}/${s}_R3_BC.sam

if [[ -f "${fastq_dir}/${s}_R1_BC_cov.fq.gz" && -f "${fastq_dir}/${s}_R3_BC_cov.fq.gz" ]]
then
        echo ${s}" has been processed."
        rm ${fastq_dir}/${s}_R1_BC.sam
        rm ${fastq_dir}/${s}_R3_BC.sam
fi

trim_galore -q 20 -j $threads --paired ${fastq_dir}/${s}_R1_BC_cov.fq.gz ${fastq_dir}/${s}_R3_BC_cov.fq.gz -o ${trim_dir}

# ### mapping
(bwa mem -SP5M -T0 -t $threads ${ref} ${trim_dir}/${s}_R1_BC_cov_val_1.fq.gz ${trim_dir}/${s}_R3_BC_cov_val_2.fq.gz | samtools view -bhS - > ${map_dir}/${s}_${genome}.bam) 2>${map_dir}/${s}_${genome}.log

### merge filesß
if [[ -f "${old_map_dir}/${s}_${genome}.bam" ]]
then
        echo "previous sequencing files found."
        samtools merge -@ $threads ${map_dir}/${s}_merged_${genome}.bam ${map_dir}/${s}_${genome}.bam ${old_map_dir}/${s}_${genome}.bam
        s=${s}_merged
fi

### identify valid pair interactions using pairtools
### all possible complex ligations are preserved (--walks-policy all)
samtools view -h ${map_dir}/${s}_${genome}.bam | \
pairtools parse --min-mapq 40 --walks-policy all --nproc-in $threads --nproc-out $threads --max-inter-align-gap 30 --chroms-path ${chrsize} --assembly ${genome} --output-stats ${map_dir}/${s}_${genome}.pairparse.txt | \
pairtools sort --nproc $threads --tmpdir=${map_dir} > ${map_dir}/${s}_${genome}.sorted.pairsam
## pairtools add col for single cell dedup. 
perl ${script_dir}/phc.pairsam_add_bc.v2.pl ${map_dir}/${s}_${genome}.sorted.pairsam
pairtools dedup --nproc-in $threads --nproc-out $threads --extra-col-pair "barcode1" "barcode2" --mark-dups --output-stats ${map_dir}/${s}_${genome}.sc.pairdedup.txt ${map_dir}/${s}_${genome}.sorted.pairsam | \
pairtools split --nproc-in $threads --nproc-out $threads --output-pairs ${map_dir}/${s}_${genome}.sc.pairs --output-sam - | samtools view -bS -@ $threads | samtools sort -T ${map_dir} -@ $threads -o ${map_dir}/${s}_${genome}.sc.pairtools.bam 

Rscript ${script_dir}/phc.summarize_pairs_lec.R ${map_dir}/${s}_${genome}.sc.pairdedup.txt
Rscript ${script_dir}/phc.plot_fragment.R ${map_dir}/${s}_${genome}.sc.pairs

if [ -f ${map_dir}/${s}_${genome}.sc.pairs ]
then
        rm ${map_dir}/${s}_${genome}.sorted.pairsam
fi

bgzip -f ${map_dir}/${s}_${genome}.sc.pairs
pairix -f ${map_dir}/${s}_${genome}.sc.pairs.gz

### count valid pair per barcode (valid reads type?)
python ${script_dir}/phc.count_pairs_sc.py --input ${map_dir}/${s}_${genome}.sc.pairs.gz --output ${map_dir}/${s}_${genome}.PairCount

### binning
bs=5000
cooler cload pairix ${chrsize}:${bs} ${map_dir}/${s}_${genome}.sc.pairs.gz ${mtx_dir}/${s}_${genome}_${bs}.cool

### Balancing
cooler zoomify --balance --balance-args '–nproc 8' -p $threads -o ${mtx_dir}/${s}_${genome}.mcool -r 5000,10000,25000,50000,100000,250000,500000,1000000,2500000 ${mtx_dir}/${s}_${genome}_${bs}.cool

