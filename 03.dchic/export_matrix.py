#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np
import glob
from collections import defaultdict
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="Generate BOTH Raw and Downsampled dcHiC matrices in a single pass.")
    parser.add_argument('--groups', nargs='+', default=['FMC', 'SMC', 'CMC'],
                        help="List of groups to process. Default: FMC SMC CMC")
    parser.add_argument('--resolutions', nargs='+', type=int, default=[10000, 50000, 100000],
                        help="Resolutions to process. Default: 10000 50000 100000")
    parser.add_argument('--raw_dir', type=str, default='../hicluster/raw_matrix',
                        help="Directory containing single-cell raw matrix files.")
    parser.add_argument('--cluster_dir', type=str, default='../cluster_chunks',
                        help="Directory containing barcode cluster chunks.")
    parser.add_argument('--chrom_sizes', type=str, default='~/genome/mm10/mm10_noYM.chrom.sizes',
                        help="Path to chrom.sizes file.")
    parser.add_argument('--out_raw', type=str, default='dchic_raw',
                        help="Output directory for full/raw matrices.")
    parser.add_argument('--out_ds', type=str, default='dchic_cell_downsampled',
                        help="Output directory for downsampled matrices.")
    parser.add_argument('--seed', type=int, default=2026,
                        help="Random seed for downsampling. Default: 2026")
    return parser.parse_args()

def load_barcodes(group_name, chunk_dir):
    barcodes = []
    pattern = os.path.join(chunk_dir, f"{group_name}_chunk*")
    files = glob.glob(pattern)
    for f in files:
        with open(f, 'r') as fh:
            lines = [line.strip().split('-')[0] for line in fh if line.strip()] 
            barcodes.extend(lines)
    return sorted(list(set(barcodes)))

def generate_bins(chrom_sizes_path, res):
    chrom_sizes_path = os.path.expanduser(chrom_sizes_path)
    chroms = pd.read_csv(chrom_sizes_path, sep='\t', names=['chrom', 'length'])
    chroms = chroms[~chroms['chrom'].str.contains('random|Un|alt|M')]
    
    bins_list = []
    global_id_counter = 1
    chrom_offset = {}
    valid_chrom_set = set(chroms['chrom'].values)
    
    for _, row in chroms.iterrows():
        chrom = row['chrom']
        length = row['length']
        chrom_offset[chrom] = global_id_counter
        n_bins = int(np.ceil(length / res))
        
        ids = np.arange(global_id_counter, global_id_counter + n_bins)
        chrom_df = pd.DataFrame({
            'chrom': chrom,
            'start': np.arange(0, n_bins * res, res),
            'end': np.arange(0, n_bins * res, res) + res,
            'id': ids
        })
        bins_list.append(chrom_df)
        global_id_counter += n_bins
        
    return pd.concat(bins_list, ignore_index=True), chrom_offset, valid_chrom_set

def process_group_dual(group, all_barcodes, ds_barcodes, raw_data_dir, chrom_offset, res, out_raw_dir, out_ds_dir, valid_chrom_set):

    print(f"  Aggregating {len(all_barcodes)} total cells (Downsampled to {len(ds_barcodes)}) for {group} at {res}...")
    
    matrix_store_raw = defaultdict(int)
    matrix_store_ds = defaultdict(int)
    ds_barcode_set = set(ds_barcodes)

    prefix_mode = None 
    for test_bc in all_barcodes[:10]:
        test_file = os.path.join(raw_data_dir, f"{test_bc}.contact.rmbkl.tsv.gz")
        if os.path.exists(test_file):
            try:
                df_head = pd.read_csv(test_file, sep='\t', header=None, nrows=1, usecols=[1])
                val = str(df_head.iloc[0,0])
                
                if val in valid_chrom_set: prefix_mode = 'none'
                elif f"chr{val}" in valid_chrom_set: prefix_mode = 'add_chr'
                elif val.startswith('chr') and val[3:] in valid_chrom_set: prefix_mode = 'remove_chr'
                
                print(f"    [Auto-Detect] File chrom col: '{val}' -> Strategy: {prefix_mode}")
                break
            except Exception as e:
                continue
    
    if prefix_mode is None:
        print("    Warning: Could not determine chrom format. Assuming standard.")
        prefix_mode = 'none'

    processed_count = 0
    reads_raw = 0
    reads_ds = 0
    
    for bc in all_barcodes:
        file_path = os.path.join(raw_data_dir, f"{bc}.contact.rmbkl.tsv.gz")
        if not os.path.exists(file_path): continue
        
        is_ds_cell = bc in ds_barcode_set
            
        try:
            chunk_iter = pd.read_csv(file_path, sep='\t', header=None, 
                                     usecols=[1, 2, 3, 4], names=['c1', 'p1', 'c2', 'p2'], 
                                     dtype={'c1': str, 'c2': str}, chunksize=100000)
            
            for chunk in chunk_iter:
                if prefix_mode == 'add_chr':
                    chunk['c1'] = 'chr' + chunk['c1']
                    chunk['c2'] = 'chr' + chunk['c2']
                elif prefix_mode == 'remove_chr':
                    chunk['c1'] = chunk['c1'].str.replace('chr', '')
                    chunk['c2'] = chunk['c2'].str.replace('chr', '')

                chunk = chunk[chunk['c1'].isin(chrom_offset) & chunk['c2'].isin(chrom_offset)]
                if chunk.empty: continue
                
                offset1 = chunk['c1'].map(chrom_offset)
                offset2 = chunk['c2'].map(chrom_offset)
                b1 = (offset1 + (chunk['p1'] // res)).values.astype(int)
                b2 = (offset2 + (chunk['p2'] // res)).values.astype(int)
                
                mask = b1 > b2
                b1[mask], b2[mask] = b2[mask], b1[mask]
                
                pairs = np.column_stack((b1, b2))
                unique_pairs, counts = np.unique(pairs, axis=0, return_counts=True)

                for (u1, u2), c in zip(unique_pairs, counts):
                    matrix_store_raw[(u1, u2)] += c
                    reads_raw += c
                    if is_ds_cell:
                        matrix_store_ds[(u1, u2)] += c
                        reads_ds += c
                        
            processed_count += 1
            if processed_count % 50 == 0:
                print(f"    Processed {processed_count}/{len(all_barcodes)} cells... (Raw Reads: {reads_raw:,} | DS Reads: {reads_ds:,})", end='\r')
                
        except Exception as e:
            pass

    print(f"\n    Finished. Total Reads - Raw: {reads_raw:,} | Downsampled: {reads_ds:,}")
    
    if reads_raw == 0:
        print("    !! ALARM: Matrix is empty! Check chrom names matching again.")
        return None, None

    mat_file_raw = os.path.join(out_raw_dir, f"{group}_{res}.matrix")
    with open(mat_file_raw, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for (b1, b2), count in matrix_store_raw.items():
            writer.writerow([b1, b2, count])

    mat_file_ds = os.path.join(out_ds_dir, f"{group}_{res}.matrix")
    with open(mat_file_ds, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for (b1, b2), count in matrix_store_ds.items():
            writer.writerow([b1, b2, count])
            
    return mat_file_raw, mat_file_ds

if __name__ == '__main__':
    args = parse_args()
    
    np.random.seed(args.seed)
    outdir_raw = os.path.abspath(args.out_raw)
    outdir_ds = os.path.abspath(args.out_ds)
    os.makedirs(outdir_raw, exist_ok=True)
    os.makedirs(outdir_ds, exist_ok=True)

    print("=== Step 1: Calculating Cell Numbers for Normalization ===")
    group_all_barcodes = {}

    for group in args.groups:
        bcs = load_barcodes(group, args.cluster_dir)
        group_all_barcodes[group] = np.array(bcs)
        print(f"  Group {group}: {len(bcs)} cells found.")

    min_cells = min(len(bcs) for bcs in group_all_barcodes.values())
    print(f"\n  [Config] Downsampling all groups to minimum size: {min_cells} cells")
    
    group_ds_barcodes = {}
    for group in args.groups:
        original_bcs = group_all_barcodes[group]
        if len(original_bcs) > min_cells:
            selected = np.random.choice(original_bcs, size=min_cells, replace=False)
        else:
            selected = original_bcs
        group_ds_barcodes[group] = selected
        print(f"  -> {group}: Selected {len(selected)} cells for DS.")

    for res in args.resolutions:
        print(f"\n=== Processing Resolution: {res} ===")
        res_dir_raw = os.path.join(outdir_raw, str(res))
        res_dir_ds = os.path.join(outdir_ds, str(res))
        os.makedirs(res_dir_raw, exist_ok=True)
        os.makedirs(res_dir_ds, exist_ok=True)
        
        bins_df, chrom_offset, valid_chrom_set = generate_bins(args.chrom_sizes, res)
        input_data_raw = []
        input_data_ds = []
        
        for group in args.groups:
            all_bcs = group_all_barcodes[group]
            ds_bcs = group_ds_barcodes[group]

            bed_file_raw = os.path.join(res_dir_raw, f"{group}_{res}.bed")
            bed_file_ds = os.path.join(res_dir_ds, f"{group}_{res}.bed")
            bins_df.to_csv(bed_file_raw, sep='\t', header=False, index=False)
            bins_df.to_csv(bed_file_ds, sep='\t', header=False, index=False)

            mat_raw, mat_ds = process_group_dual(group, all_bcs, ds_bcs, args.raw_dir, chrom_offset, res, res_dir_raw, res_dir_ds, valid_chrom_set)
            
            if mat_raw and mat_ds:
                input_data_raw.append({
                    'mat': os.path.abspath(mat_raw), 
                    'bed': os.path.abspath(bed_file_raw), 
                    'rep': f"{group}_{res}", 
                    'exp': group
                })
                input_data_ds.append({
                    'mat': os.path.abspath(mat_ds), 
                    'bed': os.path.abspath(bed_file_ds), 
                    'rep': f"{group}_{res}_ds", 
                    'exp': group
                })
            
        if input_data_raw:
            df_raw = pd.DataFrame(input_data_raw)
            txt_raw = os.path.join(outdir_raw, f"input_full_{res}.txt")
            df_raw.to_csv(txt_raw, sep='\t', header=False, index=False)
            print(f"Generated Raw Input: {txt_raw}")
            
        if input_data_ds:
            df_ds = pd.DataFrame(input_data_ds)
            txt_ds = os.path.join(outdir_ds, f"input_cell_downsampled_{res}.txt")
            df_ds.to_csv(txt_ds, sep='\t', header=False, index=False)
            print(f"Generated DS Input: {txt_ds}")