#!/bin/bash
#PBS -q hotel
#PBS -N schic_gs_from_cool
#PBS -l nodes=1:ppn=16
#PBS -l walltime=48:00:00

if [ -z ${r+x} ]; then echo "resolution for gene score calculation is not specified. Default to 10kb"; r=10; fi

#######
current=/projects/ps-renlab2/y2xie/projects/77.LC/62.schiC_MouseCortex_NovaSeq_230723/
imputed_matrix=${current}/hicluster/imputed_matrix
cell_table=${imputed_matrix}/${r}kb_resolution/filelist/cell_table.tsv ### two columns: cell id, cool url
#######

if [ -z ${genome+x} ]; then echo "genome is not specified. Default to mm10"; genome="mm10"; fi
if [[ "${genome}" == "mm10" ]]
then
    gene_meta=/projects/ps-renlab/y2xie/projects/genome_ref/Paired-Tag/mm10/mm10.gcode.vm25.bed
    chrom_size=/projects/ps-renlab2/y2xie/projects/genome_ref/mm10.main.chrom.sizes
elif [[ "${genome}" == "hg38" ]]
then
    gene_meta=/projects/ps-renlab/y2xie/projects/genome_ref/Paired-Tag/hg38/hg38.gcode.10X.bed
    chrom_size=/projects/ps-renlab2/y2xie/projects/genome_ref/hg38.main.chrom.sizes
fi

if [[ -f ${cell_table} ]]
then
	echo "calculate gene score..."
else
	ls ${imputed_matrix}/${r}kb_resolution/cool/ | grep '.cool' | awk -v p="${imputed_matrix}/${r}kb_resolution/cool" '{printf("%s/%s\n", p, $0)}' > ${imputed_matrix}/${r}kb_resolution/filelist/cell.txt
	paste <(awk -F'/' '{print $NF}' ${imputed_matrix}/${r}kb_resolution/filelist/cell.txt | cut -d. -f1) ${imputed_matrix}/${r}kb_resolution/filelist/cell.txt | sort -k1,1 > ${imputed_matrix}/${r}kb_resolution/filelist/cell.txt.tmp 
	mv ${imputed_matrix}/${r}kb_resolution/filelist/cell.txt.tmp ${imputed_matrix}/${r}kb_resolution/filelist/cell.txt
fi

mkdir ${imputed_matrix}/${r}kb_resolution/genescore
/projects/ps-renlab/y2xie/anaconda3/envs/schicluster/bin/hicluster gene-score --cell_table ${cell_table} --gene_meta ${gene_meta} --resolution ${r}000 --output_hdf_path ${imputed_matrix}/${r}kb_resolution/genescore/geneimputescore.hdf --chrom_size_path ${chrom_size} --cpu 16 --mode impute
