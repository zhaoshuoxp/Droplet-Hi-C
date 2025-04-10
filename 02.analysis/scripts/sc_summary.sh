#!/bin/bash

usage() {
cat <<EOF
Usage: sc_summary.sh [-h] [-i input.mcool] [-g genome] [-b binsize] [-o out_path]
Description:
Options:
    -h           Print help and exit
    -i           input mcool file
    -g           genome used. Support hg38 / mm10. 
    -b           binsize, default as 1Mb.
    -o           output path. Will be created if not exists
EOF
    exit 1
}

while getopts ":i:g:b:o:h" flag
do
    case "${flag}" in
        i) f=${OPTARG};;
        g) genome=${OPTARG};;
		b) bin=${OPTARG};;
        o) out_path=${OPTARG};;
        h) usage;;
    esac
done

if [ -z "${f}" ] || [ -z "${genome}" ] || [ -z "${out_path}" ]; then
    usage
fi

if [ $genome == "mm10" ]
then
	chrsize=/home/quanyiz/genome/mm10/mm10.chrom.sizes
	species=mouse
elif [ $genome == "hg38" ] 
then
	chrsize=/home/quanyiz/genome/hg38/hg38.chrom.sizes
	species=human
fi

if [ -z "${bin+x}" ]; then
	bin=1000000
fi

fname=`basename $f .mcool`
indir=`dirname $f`
mkdir -p ${out_path}/${fname} ${out_path}/tmp

/projects/ps-renlab/y2xie/anaconda3/envs/EagleC/bin/calculate-cnv -H ${f}::resolutions/${bin} -g ${genome} -e Arima --output ${out_path}/${fname}/${bin}.CNV.bedGraph --cachefolder ${out_path}/tmp
/projects/ps-renlab/y2xie/anaconda3/envs/EagleC/bin/plot-cnv --cnv-profile ${out_path}/${fname}/${bin}.CNV.bedGraph --output-figure-name ${out_path}/${fname}/${bin}.CNV.png --dot-size 0.5 --dot-alpha 0.2 --line-width 1 --boundary-width 0.5 --label-size 6 --tick-label-size 4 # --maximum-value 5 --minimum-value -5 --clean-mode

/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin/Rscript /projects/ps-renlab2/y2xie/scripts/Paired-HiC/neoloop_filter_residual.R ${out_path}/${fname}/${bin}.CNV.bedGraph ${out_path}/${fname}
/projects/ps-renlab/y2xie/anaconda3/envs/EagleC/bin/cooler dump -t pixels -H --fill-lower --join --na-rep 0 ${f}::/resolutions/${bin} > ${out_path}/${fname}/matrix.mtx
/projects/ps-renlab/y2xie/anaconda3/envs/EagleC/bin/cooler dump -t bins --na-rep 0 -H ${f}::/resolutions/${bin} > ${out_path}/${fname}/barcode.tsv
/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin/Rscript /projects/ps-renlab2/y2xie/scripts/Paired-HiC/scMatrix_summary.R ${out_path}/${fname} ${out_path}/${fname}/${bin}.CNV.bedGraph ${out_path}/${fname}
