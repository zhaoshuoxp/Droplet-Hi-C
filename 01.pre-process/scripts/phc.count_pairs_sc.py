#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='filter bam based on barcode name')
parser.add_argument('--input', type=str, dest="iput", help='.pairs or .pairs.gz input')
parser.add_argument('--output', type=str, dest="oput", help='output .stat.csv prefix')

args = parser.parse_args()

import os
import numpy as np
import pandas as pd
from time import perf_counter as pc
import gzip

def run():
    """ Run standard NMF on rank """
    start_time = pc()
    """ init input files """
    iiput = str(args.iput)
    ooput = str(args.oput)
    print("count contacts per cell...")
    count_contacts(iiput, ooput)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def count_contacts(iiput, ooput):
    cid = {}
    with oppf(iiput, 'r') as infile:
        for line in infile:
            try:
                dline = line.decode('utf8')
            except AttributeError:
                dline = line
                pass
            if dline[:1] == "#":
                continue
            readname, chrom1, pos1, chrom2, pos2, strand1, strand2, pair_type, cell1, cell2  = dline.strip("\n").split("\t")[:10]
            if cell1 in cid:
                cid[cell1]['total'] += 1
            else: ### initialize
                cid[cell1] = {'total': 1, 'mapped': 0, 'unmapped': 0, 'duplicate': 0, 'cis': 0, 'cis_1kb-': 0, 'cis_1kb+': 0, 'cis_10kb+': 0, 'trans': 0}
            if pair_type.upper() in ["NN","XX"]: 
                cid[cell1]['unmapped'] += 1
            elif chrom1 != "!" and chrom2 != "!": ### should be equal with above?
                cid[cell1]['mapped'] += 1
                if pair_type.upper() == "DD":
                    cid[cell1]['duplicates'] += 1
                elif chrom1 == chrom2:
                    cid[cell1]['cis'] += 1
                    dist = np.abs(int(pos2) - int(pos1))
                    if dist < 1000:
                        cid[cell1]['cis_1kb-'] += 1
                    if dist > 1000:
                        cid[cell1]['cis_1kb+'] += 1
                    if dist > 10000:
                        cid[cell1]['cis_10kb+'] += 1
                else:
                    cid[cell1]['trans'] += 1
            
    cidd = pd.DataFrame.from_dict(cid)
    ### filtering beofre save
    # cidd = cidd.T.loc[cidd['total'] >= 50]
    cidd.T.to_csv(ooput+".stat.csv", sep='\t')

def oppf(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode) 
    else:
        return open(filename, mode)

if __name__ == "__main__":
    run()


