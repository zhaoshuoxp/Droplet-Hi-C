#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='extracts pairs.gz based on barcode name')
parser.add_argument('--indir', type=str, dest="iput", help='input directory')
parser.add_argument('--cluster', type=str, dest="cluster", help='four columns meta: sample.pairs.gz, library, barcode, cluster')
parser.add_argument('--max', type=int, dest="maxf", default=1000, help='Max open files at the same time?')
parser.add_argument('--outdir', type=str, dest="oput", help='output directory')

args = parser.parse_args()

import os
import numpy as np
import pandas as pd
from time import perf_counter as pc
import gzip
import resource as res

def run():
    dirPrefix = str(args.iput)
    clusterf = args.cluster
    outPrefix = str(args.oput)
    max_open_files = int(args.maxf)
    start_time = pc()
    print("split contacts...")
    split_contacts(dirPrefix, outPrefix, clusterf, max_open_files)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def split_contacts(dirPrefix, outPrefix, clusterf, max_open_files):
    clust_dat = pd.read_csv(clusterf, sep="\t")
    clust_dat["uniq_barcode"] = clust_dat.apply(lambda row: f"{row['library']}:{row['barcode']}", axis=1)
    clust_dat_dict = pd.Series(clust_dat["cluster"].values, index=clust_dat["uniq_barcode"]).to_dict()
    outf_dict = dict()
    open_files = {}

    for cls_id in pd.unique(clust_dat["cluster"]):
        outf_dict[cls_id] = "".join((outPrefix, "/", str(cls_id), ".pairs"))
    
    ### writer header at once to prevent IO overhead
    print("start writing header...")
    sample_id = pd.unique(clust_dat["sample"])
    samplef = "/".join((dirPrefix, sample_id[0]))
    header_lines = []
    # Read and save header
    with oppf(samplef, 'rt') as infile:
        for dline in infile:
            if dline.startswith("#"):
                header_lines.append(dline)
            else:
                break
    header = "".join(header_lines)
    for ofile in outf_dict.values():
        with open(ofile, "a") as p:
            p.write(header)
    print("finish writing header...")

    ### split cells by cluster
    for sample in sample_id:
        library_id = pd.unique(clust_dat.loc[clust_dat['sample'] == sample]['library'])
        library = library_id[0]
        split_contacts_worker(clust_dat_dict, outf_dict, sample, library, dirPrefix, outPrefix, open_files, max_open_files)

def split_contacts_worker(clust_dat_dict, outf_dict, sample, library, dirPrefix, outPrefix, open_files, max_open_files):
    samplef = "/".join((dirPrefix, sample))
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
                wline = '\t'.join((fields))
                cls = clust_dat_dict[cid]
                if cls not in open_files:
                    open_files[cls] = open(outf_dict[cls], "a")
                open_files[cls].write(wline + "\n")
                # ofile = outf_dict[cls]
                # with open(ofile, "a") as p:
                #     p.write(wline + "\n")
            if len(open_files) >= max_open_files:
                for file in open_files.values():
                    file.close()
                open_files.clear()

def oppf(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode) 
    else:
        return open(filename, mode)

if __name__ == "__main__":
    run()


