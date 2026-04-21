#!/usr/bin/env python3
import os
import pandas as pd
import glob

# ================= Configuration =================
input_dir = 'dchic_raw'  # Directory containing input files
bias_dir = 'biases'      # Output directory for bias files
os.makedirs(bias_dir, exist_ok=True)

# ================= Main Logic =================

# 1. Automatically find all input files for different resolutions
input_files = glob.glob(os.path.join(input_dir, 'input_full_*.txt'))

if not input_files:
    print(f"Warning: No input_full_*.txt files found in {input_dir}! Please check the path.")
else:
    print(f"Found input files: {[os.path.basename(f) for f in input_files]}")
    print("Starting batch generation of Dummy Biases (All 1.0)...\n")

# 2. Outer loop: Iterate through each input file (resolution)
for raw_input_file in input_files:
    print(f"====== Processing Input File: {raw_input_file} ======")
    
    try:
        # Read input_full_*.txt (No header)
        # Format: matrix_path, bed_path, sample_name, group_name
        df = pd.read_csv(raw_input_file, sep='\t', header=None)
        
        # 3. Inner loop: Iterate through each sample in the file
        for _, row in df.iterrows():
            try:
                bed_path = row[1]       # Column 2 is BED path
                sample_name = row[2]    # Column 3 is sample name (e.g., FMC_50000)
                
                # Read BED file (chr, start, end, id)
                bed_df = pd.read_csv(bed_path, sep='\t', header=None)
                
                # Generate Bias DataFrame
                # FitHiC Format: chr, mid, bias
                dummy_bias = pd.DataFrame()
                dummy_bias[0] = bed_df[0]                     # chr
                dummy_bias[1] = (bed_df[1] + bed_df[2]) // 2  # mid point
                dummy_bias[2] = 1.0                           # <--- Force to 1.0
                
                # Save
                out_file = os.path.join(bias_dir, f"{sample_name}.biases")
                dummy_bias.to_csv(out_file, sep='\t', header=False, index=False)
                
                # Gzip
                os.system(f"gzip -f {out_file}")
                print(f"  [OK] Saved: {out_file}.gz")
                
            except Exception as e:
                print(f"  [Error] Failed processing sample {sample_name}: {e}")

    except Exception as e:
        print(f"  [Error] Failed reading input file {raw_input_file}: {e}")

print("\nDone. Biases reset to 1.0 for all resolutions.")