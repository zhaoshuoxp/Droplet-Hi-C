#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
FITHIC_PATH=$(which fithic)
PYTHON_PATH=$(which python3)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "=========================================="
echo " STEP 1: export matrix"
echo "=========================================="
python "$SCRIPT_DIR/export_matrix.py" --resolutions 10000 50000 100000

echo "Calculating biases..."
python "$SCRIPT_DIR/dumpy_bias.py" 

cp biases/FMC_10000.biases.gz biases/FMC_10000_ds.biases.gz
cp biases/SMC_10000.biases.gz biases/SMC_10000_ds.biases.gz
cp biases/CMC_10000.biases.gz biases/CMC_10000_ds.biases.gz

echo "=========================================="
echo " STEP 2: (Raw) - Compartment "
echo "=========================================="
Rscript ~/baldar/app/dcHiC/dchicf.r \
    --file input_10000.txt \
    --pcatype analyze \
    --dirovwt T \
    --diffdir Result_10000 \
    --genome mm10

echo "=========================================="
echo " STEP 3:  (Raw) - Loop "
echo "=========================================="
mkdir -p DifferentialResult/Result_10000_Loop

Rscript ~/baldar/app/dcHiC/dchicf.r \
    --file dchic_raw/input_full_10000.txt \
    --pcatype fithic \
    --diffdir Result_10000_Loop \
    --fithicpath "$FITHIC_PATH" \
    --pythonpath "$PYTHON_PATH" \
    --dirovwt T \
    --fdr 0.1

mkdir -p DifferentialResult/Result_10000
ln -s -f $(pwd)/DifferentialResult/Result_10000_Loop/fithic_run/ DifferentialResult/Result_10000/

Rscript ~/baldar/app/dcHiC/dchicf.r \
    --file dchic_raw/input_full_10000.txt \
    --pcatype cis \
    --dirovwt T 

Rscript ~/baldar/app/dcHiC/dchicf.r \
    --file dchic_raw/input_full_10000.txt \
    --pcatype dloop \
    --diffdir Result_10000 \
    --dirovwt T \
    --maxd 2000000 \
    --fdr 0.1

echo "=========================================="
echo " STEP 4:  (Downsampled) - Loop "
echo "=========================================="

mkdir -p DifferentialResult/Result_ds_10000_Loop
Rscript ~/baldar/app/dcHiC/dchicf.r \
    --file dchic_cell_downsampled/input_cell_downsampled_10000.txt \
    --pcatype fithic \
    --diffdir Result_ds_10000_Loop \
    --fithicpath "$FITHIC_PATH" \
    --pythonpath "$PYTHON_PATH" \
    --dirovwt T \
    --fdr 0.1

echo "Making edgeR matrix..."
python "$SCRIPT_DIR/build_loop_matrix.py" --samples SMC FMC CMC --res 10000

echo "=========================================="
echo " Pipeline done！"
echo " - Compartment : DifferentialResult/Result_10000_Compartment"
echo " - Raw Loop : DifferentialResult/Result_10000"
echo " - Downsampled Loop : DifferentialResult/Result_ds_10000"
echo " - edgeR Matrix : DifferentialResult/edgeR"
echo "=========================================="