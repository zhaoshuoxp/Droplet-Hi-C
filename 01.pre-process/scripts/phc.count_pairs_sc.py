#!/usr/bin/env python3

import argparse
import gzip
import pandas as pd
from time import perf_counter as pc

def parse_args():
    parser = argparse.ArgumentParser(description='Count valid pairs per barcode for scHi-C')
    parser.add_argument('--input', type=str, dest="iput", required=True, help='.pairs or .pairs.gz input')
    parser.add_argument('--output', type=str, dest="oput", required=True, help='output .stat.csv prefix')
    return parser.parse_args()

def run():
    args = parse_args()
    start_time = pc()

    print(f"Counting contacts per cell from {args.iput}...")
    count_contacts(args.iput, args.oput)

    print(f"Used (secs): {pc() - start_time:.2f}")

def count_contacts(iiput, ooput):
    cid = {}
    cols = ['total', 'mapped', 'unmapped', 'duplicate', 'cis', 'cis_1kb-', 'cis_1kb+', 'cis_10kb+', 'trans']

    open_func = gzip.open if iiput.endswith('.gz') else open

    with open_func(iiput, 'rt', encoding='utf-8') as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            sp = line.split("\t", 10)
            if len(sp) < 9:
                continue

            chrom1, pos1, chrom2, pos2 = sp[1], sp[2], sp[3], sp[4]
            pair_type = sp[7].upper()
            cell1 = sp[8]

            if cell1 not in cid:
                cid[cell1] = {
                    'total': 0, 'mapped': 0, 'unmapped': 0, 'duplicate': 0,
                    'cis': 0, 'cis_1kb-': 0, 'cis_1kb+': 0, 'cis_10kb+': 0, 'trans': 0
                }

            cell_data = cid[cell1]
            cell_data['total'] += 1

            if pair_type in ("NN", "XX"):
                cell_data['unmapped'] += 1
            elif chrom1 != "!" and chrom2 != "!":
                cell_data['mapped'] += 1

                if pair_type == "DD":
                    cell_data['duplicate'] += 1
                elif chrom1 == chrom2:
                    cell_data['cis'] += 1
                    dist = abs(int(pos2) - int(pos1))

                    if dist < 1000:
                        cell_data['cis_1kb-'] += 1
                    else:
                        cell_data['cis_1kb+'] += 1

                    if dist > 10000:
                        cell_data['cis_10kb+'] += 1
                else:
                    cell_data['trans'] += 1

    print("Formatting and saving data...")
    cidd = pd.DataFrame.from_dict(cid, orient='index').reindex(columns=cols)
    cidd.index.name = 'barcode'

    cidd.to_csv(ooput + ".stat.csv", sep='\t')

if __name__ == "__main__":
    run()
