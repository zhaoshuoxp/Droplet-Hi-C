#!/bin/bash

#######
current=$(pwd)
aln_dir=${current}/mapping
map_dir=${current}/splitPairs
in_dir=${current}/hicluster
raw_matrix=${in_dir}/raw_matrix
cell_matrix=${in_dir}/cell_matrix
imputed_matrix=${in_dir}/imputed_matrix

genome="mm10"
threads=24
res_array=(10)

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -s|--sample) s="$2"; shift ;;
        -t|--threads) threads="$2"; shift ;;
        -g|--genome) genome="$2"; shift ;;
        -r|--res) IFS=',' read -ra res_array <<< "$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$s" ]]; then
    echo "Sample prefix (-s) is required!"
    exit 1
fi

#######

mkdir -p ${map_dir} ${cell_matrix} ${raw_matrix} ${imputed_matrix}

if [[ "${genome}" == "mm10" ]]
then
    chrom_size="/home/quanyiz/genome/mm10/mm10.chrom.sizes"
    bl="/home/quanyiz/genome/mm10/mm10-blacklist.v2.bed"
elif [[ "${genome}" == "hg38" ]]
then
    chrom_size="/home/quanyiz/genome/hg38/hg38.chrom.sizes"
    bl="/home/quanyiz/genome/hg38/hg38-blacklist.v2.bed" 
fi

cd ${current}

echo -e "barcode\tlibrary\tsample\tcluster" > cluster.meta
awk -v SAMPLE=${s} '{print $1"\t"SAMPLE"\t"SAMPLE"_mm10.sc.pairs.gz\t"$1}' ${s}_cutoff.cells.txt >> cluster.meta
ulimit -n 10000
python /nfs/baldar/quanyiz/app/Droplet-Hi-C/01.pre-process/scripts/phc.sc_splitPairs_v3.py --indir ${aln_dir} --cluster cluster.meta --outdir ${map_dir} --threads $threads --max 2048

ls -1 ${map_dir}/*.pairs  > contact_table.txt
paste <(awk -F'/' '{print $NF}' contact_table.txt | cut -d. -f1) contact_table.txt | sort -k1,1 > contact_table.tsv

cell_table=${current}/contact_table.tsv
start=1
end=$(wc -l $cell_table|awk '{print $1}')

process_cell() {
    i=$1
    s=$(awk -v iline=${i} 'NR==iline {print $1}' ${cell_table})
    echo "[`date`] Start cell: $s (line $i)"

    awk -v iline=${i} 'NR==iline {print $0}' ${cell_table} > ${raw_matrix}/${s}.txt

    hicluster filter-contact --output_dir ${raw_matrix} --blacklist_1d_path ${bl} \
        --chr1 1 --pos1 2 --chr2 3 --pos2 4 \
        --contact_table ${raw_matrix}/${s}.txt \
        --chrom_size_path ${chrom_size} \
        --not_remove_duplicates

    for r in "${res_array[@]}"; do
        mkdir -p ${cell_matrix}/${r}kb_resolution ${imputed_matrix}/${r}kb_resolution/cool

        hicluster generatematrix-cell --infile ${raw_matrix}/${s}.contact.rmbkl.tsv.gz \
            --outdir ${cell_matrix}/${r}kb_resolution/ \
            --chrom_file ${chrom_size} --res ${r}000 --cell ${s} \
            --chr1 1 --pos1 2 --chr2 3 --pos2 4

        cd ${cell_matrix}/${r}kb_resolution
        for c in chr*; do
            [[ "$c" == "chr*" ]] && break
            [[ "$c" == "chrM" || "$c" == "chrY" ]] && continue

            mkdir -p "${imputed_matrix}/${r}kb_resolution/${c}"

            pad=2
            [[ "$r" == "100" ]] && pad=1

            hicluster impute-cell --indir "${cell_matrix}/${r}kb_resolution/${c}/" \
                --outdir "${imputed_matrix}/${r}kb_resolution/${c}/" \
                --cell "$s" --chrom "$c" --res "${r}000" --chrom_file "$chrom_size" \
                --output_dist 500000000 --window_size 500000000 --step_size 500000000 \
                --pad "$pad" --std 1 --rp 0.5 --tol 0.01

            out_dir="${imputed_matrix}/${r}kb_resolution/${c}/"
            for output_file in "$out_dir"*; do
                [[ "$output_file" == *.hdf5 ]] && ln -sf "$output_file" "${output_file%?}"
            done
        done
        cd ${current}

        hic-internal aggregate-chromosomes \
            --chrom_size_path ${chrom_size} \
            --resolution ${r}000 \
            --input_dir ${imputed_matrix}/${r}kb_resolution/ \
            --output_path ${imputed_matrix}/${r}kb_resolution/cool/${s}.cool \
            --chrom_wildcard "{chrom}/${s}_{chrom}_pad${pad}_std1.0_rp0.5_sqrtvc.hdf" \
            --csr
    done

    echo "[`date`] Finish cell: $s"
}

run_pool() {
    while (( $(jobs -rp | wc -l) >= threads )); do
        sleep 1
    done
    process_cell "$@" &
}

for ((i=start; i<end; i++)); do
    run_pool "$i"
done

wait
echo "✅ All jobs finished."
