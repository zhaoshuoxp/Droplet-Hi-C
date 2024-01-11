#!/bin/bash
#PBS -q home
#PBS -N schic_prepare_imputation_${start}_${end}
#PBS -l nodes=1:ppn=8
#PBS -l walltime=124:00:00

### 1000 jobs per batch

#######
current=/projects/ps-renlab2/y2xie/projects/77.LC/62.schiC_MouseCortex_NovaSeq_230723
map_dir=${current}/03.mapping/single_cell
in_dir=${current}/hicluster
raw_matrix=${in_dir}/raw_matrix
cell_matrix=${in_dir}/cell_matrix
imputed_matrix=${in_dir}/imputed_matrix
cell_table=${current}/hicluster/cell_table.txt
#######

mkdir -p ${cell_matrix} ${raw_matrix} ${imputed_matrix}

if [ -z ${genome+x} ]; then echo "genome is not specified. Default to mm10"; genome="mm10"; fi
if [[ "${genome}" == "mm10" ]]
then
    chrom_size=/projects/ps-renlab2/y2xie/projects/genome_ref/mm10.noXY.chrom.sizes
    bl=/projects/ps-renlab/y2xie/projects/genome_ref/mm10-blacklist.v2.bed
elif [[ "${genome}" == "hg38" ]]
then
    chrom_size=/projects/ps-renlab2/y2xie/projects/genome_ref/hg38.noXY.chrom.sizes
    bl=/projects/ps-renlab/y2xie/projects/genome_ref/hg38-blacklist.v2.bed
fi

cd ${current}

for ((i=${start}; i<${end}; i++)) ### manual batch
do
	s=$(awk -v iline=${i} 'NR==iline {print $1}' ${cell_table})
	
	### filter blacklist
	/projects/ps-renlab/y2xie/anaconda3/envs/schicluster/bin/hicluster filter-contact --output_dir ${raw_matrix} --blacklist_1d_path ${bl} --chr1 1 --pos1 2 --chr2 3 --pos2 4 --contact_table ${raw_matrix}/${s}.txt --chrom_size_path ${chrom_size} --not_remove_duplicates
	for r in 10 #25 100
	do
		mkdir -p ${cell_matrix}/${r}kb_resolution ${imputed_matrix}/${r}kb_resolution/cool
		egrep -v '#' ${map_dir}/${s}.pairs > ${raw_matrix}/${s}.txt
		/projects/ps-renlab/y2xie/anaconda3/envs/schicluster/bin/hicluster generatematrix-cell --infile ${raw_matrix}/${s}.txt --outdir ${cell_matrix}/${r}kb_resolution/ --chrom_file ${chrom_size} --res ${r}000 --cell ${s} --chr1 1 --pos1 2 --chr2 3 --pos2 4
		cd ${cell_matrix}/${r}kb_resolution 
		for c in chr*
		do
		    if [ "$c" == "chrX" ] || [ "$c" == "chrY" ] ### euchromsomes only
		    then
		        continue
		    else
		                mkdir -p ${imputed_matrix}/${r}kb_resolution/${c}
		                if [ "$r" == "100" ]
		                then
		                         pad=1
		                        /projects/ps-renlab/y2xie/anaconda3/envs/schicluster/bin/hicluster impute-cell --indir ${cell_matrix}/${r}kb_resolution/${c}/ --outdir ${imputed_matrix}/${r}kb_resolution/${c}/ --cell ${s} --chrom ${c} --res ${r}000 --chrom_file ${chrom_size} --output_dist 500000000 --window_size 500000000 --step_size 500000000 --pad ${pad} --std 1 --rp 0.5 --tol 0.01
		                elif [ "$r" == "25" ]
		                then
		                         pad=2
		                        /projects/ps-renlab/y2xie/anaconda3/envs/schicluster/bin/hicluster impute-cell --indir ${cell_matrix}/${r}kb_resolution/${c}/ --outdir ${imputed_matrix}/${r}kb_resolution/${c}/ --cell ${s} --chrom ${c} --res ${r}000 --chrom_file ${chrom_size} --output_dist 10050000 --window_size 500000000 --step_size 500000000 --pad ${pad} --std 1 --rp 0.5 --tol 0.01
		                elif [ "$r" == "10" ]
		                then
		                        pad=2
		                        /projects/ps-renlab/y2xie/anaconda3/envs/schicluster/bin/hicluster impute-cell --indir ${cell_matrix}/${r}kb_resolution/${c}/ --outdir ${imputed_matrix}/${r}kb_resolution/${c}/ --cell ${s} --chrom ${c} --res ${r}000 --chrom_file ${chrom_size} --output_dist 5050000 --window_size 30000000 --step_size 10000000 --pad ${pad} --std 1 --rp 0.5 --tol 0.01

		                fi
		    fi
		done
		/projects/ps-renlab/y2xie/anaconda3/envs/schicluster/bin/hic-internal aggregate-chromosomes --chrom_size_path ${chrom_size} --resolution ${r}000 --input_dir ${imputed_matrix}/${r}kb_resolution/ --output_path ${imputed_matrix}/${r}kb_resolution/cool/${s}.cool --chrom_wildcard "{chrom}/${s}_{chrom}_pad${pad}_std1.0_rp0.5_sqrtvc.hdf5" --csr
	done
done
