#!/usr/bin/env python
import argparse
import h5py
import cooler
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import rankdata
from scipy.io.mmio import MMFile

def parse_args():
    parser = argparse.ArgumentParser(description="Process gene data and generate mtx for gene score.")
    parser.add_argument('-r', '--res', type=int, required=True, help="Resolution for binning.")
    parser.add_argument('--genome', type=str, required=True, help="Genome identifier.")
    parser.add_argument('chrom_size_path', type=str, help="Path to chromosome size file.")
    parser.add_argument('gene_hdf_path', type=str, help="Path to the gene HDF file.")
    parser.add_argument('sample', type=str, help="Sample name.")
    return parser.parse_args()
    
# Override MMFile to prevent scientific notation in output
class MMFileFixedFormat(MMFile):
    def _field_template(self, field, precision):
        return f'%.{precision}f\n'

def main():
    # Parse arguments
    args = parse_args()

    # Load chromosome sizes
    chrom_sizes = cooler.read_chromsizes(args.chrom_size_path, all_names=True)

    # Load gene HDF data
    gene_hdf = pd.read_hdf(args.gene_hdf_path, key='data')

    # Load gene meta data
    gene_meta_path = f'/home/quanyiz/baldar/genome/{args.genome}/genes_pri.bed'
    gene_meta = pd.read_csv(gene_meta_path, names=['chrom', 'start', 'end', 'gene_id', 'gene_name'], index_col='gene_id', sep='\t')
    gene_meta['bin_len'] = (gene_meta['end'] // args.res) - (gene_meta['start'] // args.res) + 1
    gene_meta['weight'] = (gene_meta['bin_len'] + 2) * (gene_meta['bin_len'] + 1) / 2

    # Read chromosome sizes and filter gene meta
    chromsize = pd.read_csv(args.chrom_size_path, sep='\t', header=None, index_col=0)
    gene_meta = gene_meta[gene_meta['chrom'].isin(chromsize.index)]

    # Load gene data into AnnData
    gene3c = anndata.AnnData(gene_hdf)

    # Load statistics
    stat = pd.read_csv(f'mapping/{args.sample}_{args.genome}.PairCount.stat.csv', sep='\t')
    stat.columns.values[0] = 'cool_cell'
    tmeta = gene3c.obs.merge(stat, left_index=True, right_on='cool_cell')
    tmeta.index = tmeta['cool_cell']
    tmeta = tmeta.drop(columns=['cool_cell'])
    gene3c = gene3c[tmeta.index, :]
    gene3c.obs = tmeta

    # Filter gene3c AnnData
    genefilter = ((gene3c.X > 0).sum(axis=0) > 10) & (gene3c.var.index.isin(gene_meta['gene_name']))
    gene3c = gene3c[:, genefilter]

    # Normalize and save
    gene3c.var.index.name = 'gene_name'
    gene3c.var_names_make_unique()
    sc.pp.normalize_total(gene3c, target_sum=np.median(gene3c.X.sum(axis=1)))
    gene3c.write(f'{args.sample}_adata_{args.res}.h5ad')

    # Identify highly variable genes
    tmp = np.std(gene3c.X, axis=0)
    tmp = rankdata(-tmp)
    gene3c.var['highly_variable'] = (tmp < 2000)
    gene3c = gene3c[:, gene3c.var['highly_variable']]

    # Create sparse matrix for output
    gene3c_mtx = csr_matrix(gene3c.X)
    gene3c_cells = gene3c.obs.index
    gene3c_genes = gene3c.var.index

    # Save barcodes and genes
    pd.DataFrame(gene3c_cells).to_csv(f'hicluster/imputed_matrix/{args.res}kb_resolution/genescore/geneimputescore.barcodes.tsv', sep='\t')
    pd.DataFrame(gene3c_genes).to_csv(f'hicluster/imputed_matrix/{args.res}kb_resolution/genescore/geneimputescore.genes.tsv', sep='\t')

    # Write the sparse matrix to a file
    MMFileFixedFormat().write(f'hicluster/imputed_matrix/{args.res}kb_resolution/genescore/geneimputescore.mtx', gene3c_mtx, precision=2)

if __name__ == "__main__":
    main()
