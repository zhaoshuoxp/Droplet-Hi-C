#!/usr/bin/env python3
import os
import cooler
import pandas as pd
import numpy as np

bulk_mcool = '../matrices/ROOT_TDT_mm10.mcool' 

samples = ['FMC', 'SMC', 'CMC'] 

resolutions = [10000, 50000, 100000] 
bias_dir = 'biases'
os.makedirs(bias_dir, exist_ok=True)

def generate_bias_from_bulk(bulk_path, sample_list, res, out_dir):
    uri = f'{bulk_path}::resolutions/{res}'
    
    try:
        c = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error loading {uri}: {e}")
        return

    bins = c.bins()[:]
    
    if 'weight' not in bins.columns:
        print(f"!! Error: {uri} missing weight ")
        print(f"   run: cooler balance {uri}")
        return

    print(f"Generating bias from bulk for resolution {res}...")

    weights = bins['weight'].fillna(0).to_numpy()
    biases = np.zeros_like(weights)
    valid_mask = weights > 0
    biases[valid_mask] = 1.0 / weights[valid_mask]
    mids = (bins['start'] + bins['end']) // 2
    
    bias_df = pd.DataFrame({
        'chr': bins['chrom'],
        'mid': mids,
        'bias': biases
    })
 
    for samp in sample_list:
        sample_prefix = f"{samp}_{res}"
        out_file = os.path.join(out_dir, f"{sample_prefix}.biases")
        
        bias_df.to_csv(out_file, sep='\t', header=False, index=False)
        os.system(f"gzip -f {out_file}")
        print(f"  Saved: {out_file}.gz")

if __name__ == '__main__':
    for res in resolutions:
        generate_bias_from_bulk(bulk_mcool, samples, res, bias_dir)
