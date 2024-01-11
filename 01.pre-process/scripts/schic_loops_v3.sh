#!/bin/bash
#PBS -q hotel
#PBS -N schic_call_loop_${s}
#PBS -l nodes=1:ppn=8
#PBS -l walltime=96:00:00

r=10
subdir="DPT"
#######
current=/projects/ps-renlab2/y2xie/projects/77.LC/62.schiC_MouseCortex_NovaSeq_230723
imputed_matrix=${current}/hicluster/imputed_matrix
meta=${imputed_matrix}/${r}kb_resolution/filelist/${subdir}/${s}_bc.txt ### one column: cell_id
#######

if [ -z ${genome+x} ]; then echo "genome is not specified. Default to mm10"; genome="mm10"; fi
if [[ "${genome}" == "mm10" ]]
then
    chrom_size=/projects/ps-renlab2/y2xie/projects/genome_ref/mm10.noXY.chrom.sizes
elif [[ "${genome}" == "hg38" ]]
then
    chrom_size=/projects/ps-renlab2/y2xie/projects/genome_ref/hg38.noXY.chrom.sizes
fi

mkdir -p ${imputed_matrix}/${r}kb_resolution/matrices/${subdir} ${imputed_matrix}/${r}kb_resolution/loop/${subdir}

# loop-bkg-cell is already calculated by schic_impute_v2.sh

for c in chr{1..19}
do
    ls ${imputed_matrix}/${r}kb_resolution/${c} | grep "pad2_std1.0_rp0.5_sqrtvc.hdf5" - | grep -f ${meta} - | sed 's/.hdf5//g' - > ${imputed_matrix}/${r}kb_resolution/filelist/${subdir}/imputematrices_${s}_${c}_1.txt

    awk -v path="${imputed_matrix}/${r}kb_resolution/${c}/" '{print path $0}' ${imputed_matrix}/${r}kb_resolution/filelist/${subdir}/imputematrices_${s}_${c}_1.txt > ${imputed_matrix}/${r}kb_resolution/filelist/${subdir}/imputematrices_${s}_${c}.txt

    /projects/ps-renlab/y2xie/anaconda3/envs/schicluster/bin/hicluster loop-sumcell-chr --cell_list ${imputed_matrix}/${r}kb_resolution/filelist/${subdir}/imputematrices_${s}_${c}.txt --outprefix ${imputed_matrix}/${r}kb_resolution/matrices/${subdir}/imputematrices_${s}_${c} --res ${r}000
done

### merge loop list by cluster
/projects/ps-renlab/y2xie/anaconda3/envs/schicluster/bin/hicluster loop-mergechr --inprefix ${imputed_matrix}/${r}kb_resolution/matrices/${subdir}/imputematrices_${s} --outprefix ${imputed_matrix}/${r}kb_resolution/loop/${subdir}/imputematrices_${s} --chrom_file ${chrom_size}
