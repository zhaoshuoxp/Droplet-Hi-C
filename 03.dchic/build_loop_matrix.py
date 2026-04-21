#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import sys
import argparse


def load_fithic_loops(filepath):
    print(f">>> Loading Loop Universe from {filepath}...")
    try:
        df = pd.read_csv(filepath, sep='\t', index_col=0)
        coords = df.index.to_series().str.split('_', expand=True)
        dist = (coords[1].astype(int) - coords[3].astype(int)).abs()
        
        valid_mask = dist > 0
        valid_loops = df.index[valid_mask]
        
        print(f"  Total Loops in file: {len(df)}")
        print(f"  Loops > 0 distance (kept): {len(valid_loops)}")

        result_df = pd.DataFrame(index=valid_loops)
        return result_df
    except Exception as e:
        print(f"  [Error] reading FithicResult: {e}")
        return None

def load_sample_counts(sample, raw_dir, res, target_index):
    print(f"\n--- Processing {sample} ---")
    bed_file = os.path.join(raw_dir, str(res), f"{sample}_{res}.bed")
    mat_file = os.path.join(raw_dir, str(res), f"{sample}_{res}.matrix")
 
    if not os.path.exists(bed_file):
        bed_file = os.path.join(raw_dir, f"{sample}_{res}.bed")
        mat_file = os.path.join(raw_dir, f"{sample}_{res}.matrix")
    
    if not os.path.exists(bed_file) or not os.path.exists(mat_file):
        print(f"  [Error] Files not found for {sample}: {bed_file} or {mat_file}")
        return pd.Series(0, index=target_index, name=sample)

    print(f"  Reading BED: {bed_file}")
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'id'])
    bed_df['coord_str'] = bed_df['chrom'] + "_" + bed_df['start'].astype(str)
    bin_map = bed_df.set_index('id')['coord_str'].to_dict()

    print(f"  Reading Matrix: {mat_file}")
    mat_df = pd.read_csv(mat_file, sep='\t', header=None, names=['b1', 'b2', 'count'])

    print("  Mapping Bin IDs to Loop IDs...")
    mat_df['s1'] = mat_df['b1'].map(bin_map)
    mat_df['s2'] = mat_df['b2'].map(bin_map)

    mat_df['LoopID'] = mat_df['s1'] + "_" + mat_df['s2']

    counts_series = mat_df.set_index('LoopID')['count']
    aligned_counts = counts_series.reindex(target_index, fill_value=0)
    
    zero_count = (aligned_counts == 0).sum()
    if zero_count > 0:
        mat_df['LoopID_Swap'] = mat_df['s2'] + "_" + mat_df['s1']
        counts_swap = mat_df.set_index('LoopID_Swap')['count']
        aligned_swap = counts_swap.reindex(target_index, fill_value=0)
        
        mask = (aligned_counts == 0) & (aligned_swap > 0)
        aligned_counts[mask] = aligned_swap[mask]

    matched_count = aligned_counts.astype(bool).sum()
    print(f"  Matched {matched_count} non-zero loops out of {len(target_index)}")
    
    return aligned_counts

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge Loop Counts for Raw and Downsampled datasets.")
    parser.add_argument("--samples", nargs='+', default=['SMC', 'FMC', 'CMC'], 
                        help="List of samples (e.g., --samples SMC FMC CMC). Default: SMC FMC CMC")
    parser.add_argument("--res", type=int, default=10000, 
                        help="Resolution of the Hi-C matrices. Default: 10000")
    
    args = parser.parse_args()
    
    SAMPLES = args.samples
    RES = args.res
    HALF_RES = int(RES / 2)

    print(f"⚙️  Configuration Loaded: RES={RES}, SAMPLES={SAMPLES}")

    RUN_CONFIGS = [
        {
            "version": "Downsampled",
            "raw_dir": "dchic_cell_downsampled", 
            "loop_file": f"DifferentialResult/Result_ds_{RES}_Loop/fithic_run/FithicResult_p.01_ct2.txt",
            "output_file": "edgeR/Loop_RawCounts_Matrix_ds.txt"
        },
        {
            "version": "Raw",
            "raw_dir": "dchic_raw", 
            "loop_file": f"DifferentialResult/Result_{RES}_Loop/fithic_run/FithicResult_p.01_ct2.txt",
            "output_file": "edgeR/Loop_RawCounts_Matrix.txt"
        }
    ]

    for config in RUN_CONFIGS:
        print(f"\n{'='*60}")
        print(f"🚀 STARTING PIPELINE: {config['version']}")
        print(f"{'='*60}")
        
        out_dir = os.path.dirname(config['output_file'])
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        final_df = load_fithic_loops(config['loop_file'])
        
        if final_df is None:
            print(f"⏭️ Skipping {config['version']} due to missing loop file.\n")
            continue 

        for samp in SAMPLES:
            counts = load_sample_counts(samp, config['raw_dir'], RES, final_df.index)
            final_df[samp] = counts
            
        final_df.to_csv(config['output_file'], sep='\t')
        print(f"\n>>> Done! {config['version']} matrix saved to: {config['output_file']}\n")
        
    print("✅ All processing finished successfully!")
