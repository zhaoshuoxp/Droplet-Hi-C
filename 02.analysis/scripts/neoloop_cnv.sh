#!/bin/bash

usage() {
cat <<EOF
Usage: neoloop_cnv.sh [-h] [-i input.mcool] [-g genome] [-n Njob_predictSV] [-o out_path] [-c calculate cnv (True / False)]
Description:
Options:
    -h           Print help and exit
    -i           input mcool file
    -g           genome used. Support hg38 / mm10. 
    -n           number of threads for predictSV
    -o           output path. Will be created if not exists
    -c           Calculate cnv or not before EagleC?
EOF
    exit 1
}

while getopts ":i:g:n:o:c:h" flag
do
    case "${flag}" in
        i) cool=${OPTARG};;
        g) genome=${OPTARG};;
        n) cpu=${OPTARG};;
        o) out_path=${OPTARG};;
		c) cnv_cal=${OPTARG};;
        h) usage;;
    esac
done

if [ -z "${cool}" ] || [ -z "${genome}" ] || [ -z "${out_path}" ]; then
    usage
fi

if [ $genome == "mm10" ]
then
	chrsize=/projects/ps-renlab/y2xie/projects/genome_ref/mm10.main.chrom.sizes
	species=mouse
elif [ $genome == "hg38" ] 
then
	chrsize=/projects/ps-renlab/y2xie/projects/genome_ref/hg38.main.chrom.sizes
	species=human
fi

if [ -z "${name+x}" ]; then
	cpu=8
fi

if [ -z "${cnv_cal+x}" ]; then
	cnv_cal=True
fi

mkdir -p ${out_path}/EagleC_SV/tmp
tmp_path=${out_path}/EagleC_SV/tmp
s=`basename $cool .mcool`

if [[ "$cnv_cal" == "True" ]]
then
	echo "start calculating cnv..."
	for bin in {5000,10000,25000,50000,100000}
	do
		echo "perform cnv calculation at $bin resolution..."
		calculate-cnv -H ${cool}::resolutions/${bin} -g ${genome} -e Arima --output ${out_path}/${s}_${bin}.CNV.bedGraph --cachefolder ${tmp_path}
		segment-cnv --cnv-file ${out_path}/${s}_${bin}.CNV.bedGraph --binsize ${bin} --ploidy 2 --cbs-pvalue 5e-2 --max-dist 1 --output ${out_path}/${s}_${bin}.CNV-seg.bedGraph
		plot-cnv --cnv-profile ${out_path}/${s}_${bin}.CNV.bedGraph --cnv-segment ${out_path}/${s}_${bin}.CNV-seg.bedGraph --output-figure-name ${out_path}/${s}_${bin}.CNV.png --dot-size 0.5 --dot-alpha 0.2 --line-width 1 --boundary-width 0.5 --label-size 6 --tick-label-size 4 
		correct-cnv -H ${cool}::resolutions/${bin} --cnv-file ${out_path}/${s}_${bin}.CNV-seg.bedGraph --nproc $cpu -f
	done
fi

echo "start predicting SV..."
### run prediction in parallel
for i in $(seq 1 $cpu); do nohup predictSV --hic-5k ${cool}::resolutions/5000 --hic-10k ${cool}::resolutions/10000 --hic-50k ${cool}::resolutions/25000 -O ${out_path}/EagleC_SV/${s} -g ${genome} --balance-type CNV --output-format full --prob-cutoff-5k 0.8 --prob-cutoff-10k 0.8 --prob-cutoff-50k 0.99999 > ${tmp_path}/nohup_${i}.log & done

