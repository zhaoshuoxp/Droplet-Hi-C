#!/bin/bash
#PBS -q hotel
#PBS -N schic_domain_from_cool_${s}
#PBS -l nodes=1:ppn=8
#PBS -l walltime=72:00:00

source /home/y2xie/anaconda3/etc/profile.d/conda.sh
conda activate seurat

if [ -z ${r+x} ]; then echo "resolution for domain calling is not specified. Default to 25kb"; r=25; fi

#######
current=/projects/ps-renlab2/y2xie/projects/77.LC/62.schiC_MouseCortex_NovaSeq_230723/
imputed_matrix=${current}/hicluster/imputed_matrix/
cell_table=${imputed_matrix}/${r}kb_resolution/domain/${s}/cell_table.csv ### two columns: cell id, cool url
# bin_dir=${imputed_matrix}/bins ### where chromsome bins are stored
#######

mkdir ${imputed_matrix}/${r}kb_resolution/domain/${s}
hicluster domain --cell_table_path ${cell_table} --output_prefix ${imputed_matrix}/${r}kb_resolution/domain/${s}/${s} --resolution ${r}000 --window_size 10 --cpu 8
