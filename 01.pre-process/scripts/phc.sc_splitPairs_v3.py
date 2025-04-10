#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='extracts pairs.gz based on barcode name')
parser.add_argument('--indir', type=str, dest="iput", help='input directory')
parser.add_argument('--cluster', type=str, dest="cluster", help='four columns meta: sample.pairs.gz, library, barcode, cluster')
parser.add_argument('--max', type=int, dest="maxf", default=1000, help='Max open files at the same time?')
parser.add_argument('--outdir', type=str, dest="oput", help='output directory')
parser.add_argument('--threads', type=int, dest="threads", default=4, help='Number of threads to use')

args = parser.parse_args()

import os
import numpy as np
import pandas as pd
from time import perf_counter as pc
import gzip
import resource as res
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

def run():
    dirPrefix = str(args.iput)
    clusterf = args.cluster
    outPrefix = str(args.oput)
    max_open_files = int(args.maxf)
    num_threads = int(args.threads)
    start_time = pc()
    print("split contacts...")
    split_contacts(dirPrefix, outPrefix, clusterf, max_open_files, num_threads)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def split_contacts(dirPrefix, outPrefix, clusterf, max_open_files, num_threads):
    clust_dat = pd.read_csv(clusterf, sep="\t")
    clust_dat["uniq_barcode"] = clust_dat.apply(lambda row: f"{row['library']}:{row['barcode']}", axis=1)
    clust_dat_dict = pd.Series(clust_dat["cluster"].values, index=clust_dat["uniq_barcode"]).to_dict()
    
    # Create output directory if it doesn't exist
    os.makedirs(outPrefix, exist_ok=True)
    
    # Get unique clusters and create output files
    unique_clusters = pd.unique(clust_dat["cluster"])
    outf_dict = {cls_id: os.path.join(outPrefix, f"{cls_id}.pairs") for cls_id in unique_clusters}
    
    # Create a lock for thread-safe file writing
    file_locks = {cls_id: Lock() for cls_id in unique_clusters}
    
    ### writer header at once to prevent IO overhead
    print("start writing header...")
    sample_id = pd.unique(clust_dat["sample"])
    samplef = os.path.join(dirPrefix, sample_id[0])
    header_lines = []
    
    # Read and save header
    with oppf(samplef, 'rt') as infile:
        for dline in infile:
            if dline.startswith("#"):
                header_lines.append(dline)
            else:
                break
    header = "".join(header_lines)
    
    # Write headers to all output files
    for ofile in outf_dict.values():
        with open(ofile, "w") as p:
            p.write(header)
    print("finish writing header...")

    ### split cells by cluster using multiple threads
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        for sample in sample_id:
            library_id = pd.unique(clust_dat.loc[clust_dat['sample'] == sample]['library'])
            library = library_id[0]
            futures.append(
                executor.submit(
                    split_contacts_worker,
                    clust_dat_dict,
                    outf_dict,
                    sample,
                    library,
                    dirPrefix,
                    file_locks
                )
            )
        
        # Wait for all tasks to complete
        for future in as_completed(futures):
            future.result()

def split_contacts_worker(clust_dat_dict, outf_dict, sample, library, dirPrefix, file_locks):
    samplef = os.path.join(dirPrefix, sample)
    open_files = {}
    
    try:
        with oppf(samplef, 'rt') as infile:
            for dline in infile:
                if dline.startswith("#"):
                    continue
                fields = dline.strip("\n").split("\t")
                if fields[-1] != fields[-2]:
                    continue
                cid = ":".join((library, fields[-1]))
                if cid in clust_dat_dict:
                    fields[-2:] = [cid, cid]
                    wline = '\t'.join(fields)
                    cls = clust_dat_dict[cid]
                    ofile = outf_dict[cls]
                    
                    # Use lock to ensure thread-safe file writing
                    with file_locks[cls]:
                        if cls not in open_files:
                            open_files[cls] = open(ofile, "a")
                        open_files[cls].write(wline + "\n")
    finally:
        # Ensure all files are closed when done
        for file in open_files.values():
            file.close()

def oppf(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode) 
    else:
        return open(filename, mode)

if __name__ == "__main__":
    run()