#!/bin/bash
#PBS -q home
#PBS -N map_hic_${s}
#PBS -l nodes=1:ppn=16
#PBS -l walltime=128:00:00

################################
current="/home/y2xie/scratch/33.Paired_HiC_NovaSeq_221206"
fastq_dir="${current}/01.rawdata/"
trim_dir="${current}/02.trimmed/"
map_dir="${current}/03.mapping/"
old_map_dir="${current}/07.prev_mapping/"
mtx_dir="${current}/04.matrices/"
script_dir="/projects/ps-renlab/y2xie/scripts/Paired-HiC/"
# bin_path="/projects/ps-renlab/y2xie/anaconda3/bin/"

mm10_bwa="/projects/ps-renlab/share/bwa_indices/mm10.fa"
hg38_bwa="/projects/ps-renlab/share/bwa_indices/hg38.fa"
mix_bwa="/projects/ps-renlab/y2xie/projects/genome_ref/GRCh38_and_mm10/fasta/genome.fa"
mm10_chrsize="/projects/ps-renlab/share/STAR_indices/mm10/chrNameLength.txt"
hg38_chrsize="/projects/ps-renlab/y2xie/projects/genome_ref/hg38.main.chrom.sizes"
mix_chrsize="/projects/ps-renlab/y2xie/projects/genome_ref/GRCh38_and_mm10/star/chrNameLength.txt"
mm10_bl="/projects/ps-renlab/y2xie/projects/genome_ref/mm10.blacklist.bed"
hg38_bl="/projects/ps-renlab/y2xie/projects/genome_ref/hg38-blacklist.v2.bed" 
atac_10X="/projects/ps-renlab/y2xie/projects/genome_ref/737K-cratac-whitelist/737K-cratac-v1"
arc_10X="/projects/ps-renlab/y2xie/projects/genome_ref/737K-crarc-whitelist/737K-arc-v1"
##################################

if [ -z ${genome+x} ]; then echo "genome is not specified. Default to mm10"; genome="mm10"; fi
if [ -z ${mode+x} ]; then echo "10X kit used is not specified. Default to arc"; mode="arc"; fi

if [ $genome == "mm10" ]; then ref=${mm10_bwa}; chrsize=${mm10_chrsize}; fi
if [ $genome == "hg38" ]; then ref=${hg38_bwa}; chrsize=${hg38_chrsize}; fi
### only for species mixing experiment!!!!
if [ $genome == "mix" ]; then ref=${mix_bwa}; chrsize=${mix_chrsize}; fi
        
if [ $mode == "arc" ]; then ref_10X=${arc_10X}; else ref_10X=${atac_10X}; fi

echo "process scHiC fastq... mode: "${mode}
/projects/ps-renlab/y2xie/packages/Paired-Tag-master/reachtools/reachtools combine_hic ${mode} ${fastq_dir}/${s}
(zcat ${fastq_dir}/${s}_R1_combined.fq.gz | bowtie -x ${ref_10X} - --nofw -m 1 -v 1 -S ${fastq_dir}/${s}_R1_BC.sam) 2>${fastq_dir}/${s}_R1.log
zcat ${fastq_dir}/${s}_R3_combined.fq.gz | bowtie -x ${ref_10X} - --nofw -m 1 -v 1 -S ${fastq_dir}/${s}_R3_BC.sam

/projects/ps-renlab/y2xie/packages/Paired-Tag-master/reachtools/reachtools convert_hic2 ${fastq_dir}/${s}_R1_BC.sam
/projects/ps-renlab/y2xie/packages/Paired-Tag-master/reachtools/reachtools convert_hic2 ${fastq_dir}/${s}_R3_BC.sam

if [[ -f "${fastq_dir}/${s}_R1_BC_cov.fq.gz" && -f "${fastq_dir}/${s}_R3_BC_cov.fq.gz" ]]
then
        echo ${s}" has been processed."
        rm ${fastq_dir}/${s}_R1_BC.sam
        rm ${fastq_dir}/${s}_R3_BC.sam
fi

trim_galore -q 20 --paired ${fastq_dir}/${s}_R1_BC_cov.fq.gz ${fastq_dir}/${s}_R3_BC_cov.fq.gz -o ${trim_dir}

# ### mapping
(bwa mem -SP5M -T0 -t16 ${ref} ${trim_dir}/${s}_R1_BC_cov_val_1.fq.gz ${trim_dir}/${s}_R3_BC_cov_val_2.fq.gz | samtools view -bhS - > ${map_dir}/${s}_${genome}.bam) 2>${map_dir}/${s}_${genome}.log

### merge files
if [[ -f "${old_map_dir}/${s}_${genome}.bam" ]]
then
        echo "previous sequencing files found."
        samtools merge -@ 16 ${map_dir}/${s}_merged_${genome}.bam ${map_dir}/${s}_${genome}.bam ${old_map_dir}/${s}_${genome}.bam
        s=${s}_merged
fi

### identify valid pair interactions using pairtools
### all possible complex ligations are preserved (--walks-policy all)
samtools view -h ${map_dir}/${s}_${genome}.bam | \
pairtools parse --min-mapq 40 --walks-policy all --nproc-in 16 --nproc-out 16 --max-inter-align-gap 30 --chroms-path ${chrsize} --assembly ${genome} --output-stats ${map_dir}/${s}_${genome}.pairparse.txt | \
pairtools sort --nproc 16 --tmpdir=${map_dir} > ${map_dir}/${s}_${genome}.sorted.pairsam
## pairtools add col for single cell dedup. 
perl ${script_dir}/phc.pairsam_add_bc.v2.pl ${map_dir}/${s}_${genome}.sorted.pairsam
pairtools dedup --nproc-in 16 --nproc-out 16 --extra-col-pair "barcode1" "barcode2" --mark-dups --output-stats ${map_dir}/${s}_${genome}.sc.pairdedup.txt ${map_dir}/${s}_${genome}.sorted.pairsam | \
pairtools split --nproc-in 16 --nproc-out 16 --output-pairs ${map_dir}/${s}_${genome}.sc.pairs --output-sam - | samtools view -bS -@ 16 | samtools sort -T ${map_dir} -@ 16 -o ${map_dir}/${s}_${genome}.sc.pairtools.bam 

/home/y2xie/anaconda3/envs/seurat/bin/Rscript ${script_dir}/phc.summarize_pairs_lec.R ${map_dir}/${s}_${genome}.sc.pairdedup.txt
/home/y2xie/anaconda3/envs/seurat/bin/Rscript ${script_dir}/phc.plot_fragment.R ${map_dir}/${s}_${genome}.sc.pairs

if [ -f ${map_dir}/${s}_${genome}.sc.pairs ]
then
        rm ${map_dir}/${s}_${genome}.sorted.pairsam
fi

bgzip -@ 16 -f ${map_dir}/${s}_${genome}.sc.pairs
pairix -f ${map_dir}/${s}_${genome}.sc.pairs.gz

### count valid pair per barcode (valid reads type?)
python ${script_dir}/phc.count_pairs_sc.py --input ${map_dir}/${s}_${genome}.sc.pairs.gz --output ${map_dir}/${s}_${genome}.PairCount

### binning
bs=5000
cooler cload pairix ${chrsize}:${bs} ${map_dir}/${s}_${genome}.sc.pairs.gz ${mtx_dir}/${s}_${genome}_${bs}.cool

### Balancing
cooler zoomify --balance --balance-args '–nproc 8' -p 16 -o ${mtx_dir}/${s}_${genome}.mcool -r 5000,10000,25000,50000,100000,250000,500000,1000000,2500000 ${mtx_dir}/${s}_${genome}_${bs}.cool
