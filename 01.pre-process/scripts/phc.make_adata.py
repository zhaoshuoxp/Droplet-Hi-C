#!/usr/bin/env python
import argparse
import os

import anndata
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, diags, issparse
from scipy.stats import rankdata
from scipy.io import mmwrite

def parse_args():
    parser = argparse.ArgumentParser(description="Process gene data and generate mtx for gene score.")
    parser.add_argument('-r', '--res', type=int, required=True, help="Resolution for binning.")
    parser.add_argument('--genome', type=str, required=True, help="Genome identifier.")
    parser.add_argument('chrom_size_path', type=str, help="Path to chromosome size file.")
    parser.add_argument('gene_hdf_path', type=str, help="Path to the gene HDF file.")
    parser.add_argument('sample', type=str, help="Sample name.")
    return parser.parse_args()


def normalize_total(matrix):
    row_sums = np.asarray(matrix.sum(axis=1)).ravel()
    target_sum = float(np.median(row_sums))
    scale = np.divide(
        target_sum,
        row_sums,
        out=np.zeros_like(row_sums, dtype=float),
        where=row_sums != 0,
    )
    if issparse(matrix):
        return diags(scale).dot(matrix)
    return np.asarray(matrix) * scale[:, None]


def column_std(matrix):
    if not issparse(matrix):
        return np.std(np.asarray(matrix), axis=0)
    means = np.asarray(matrix.mean(axis=0)).ravel()
    means_squared = np.asarray(matrix.power(2).mean(axis=0)).ravel()
    return np.sqrt(np.maximum(means_squared - means**2, 0))


def main():
    # Parse arguments
    args = parse_args()

    # Load gene HDF data
    gene_hdf = pd.read_hdf(args.gene_hdf_path, key='data')

    # Load gene meta data
    gene_meta_path = os.path.join(os.path.dirname(args.chrom_size_path), 'genes_pri.bed')
    if not os.path.isfile(gene_meta_path):
        raise FileNotFoundError(f"Gene metadata is not readable: {gene_meta_path}")
    gene_meta = pd.read_csv(gene_meta_path, names=['chrom', 'start', 'end', 'gene_id', 'gene_name'], index_col='gene_id', sep='\t')

    # Read chromosome sizes and filter gene meta
    chromsize = pd.read_csv(args.chrom_size_path, sep='\t', header=None, index_col=0)
    gene_meta = gene_meta[gene_meta['chrom'].isin(chromsize.index)]

    # Load gene data into AnnData
    gene3c = anndata.AnnData(gene_hdf)

    # Load statistics
    stat = pd.read_csv(f'mapping/{args.sample}_{args.genome}.PairCount.stat.csv', sep='\t')
    stat.columns.values[0] = 'cool_cell'
    if stat['cool_cell'].duplicated().any():
        raise ValueError("Pair-count statistics contain duplicate cell identifiers")
    stat = stat.set_index('cool_cell')
    common_cells = gene3c.obs_names[gene3c.obs_names.isin(stat.index)]
    if len(common_cells) == 0:
        raise ValueError("Gene scores and pair-count statistics have no cells in common")
    tmeta = stat.loc[common_cells].copy()
    gene3c = gene3c[common_cells, :].copy()
    gene3c.obs = tmeta

    # Filter gene3c AnnData
    expressed_cells = np.asarray((gene3c.X > 0).sum(axis=0)).ravel()
    genefilter = (expressed_cells > 10) & gene3c.var_names.isin(gene_meta['gene_name'])
    gene3c = gene3c[:, genefilter].copy()

    # Normalize and save
    gene3c.var.index.name = 'gene_name'
    gene3c.var_names_make_unique()
    gene3c.X = normalize_total(gene3c.X)
    gene3c.write(f'{args.sample}_adata_{args.res}.h5ad')

    # Identify highly variable genes
    tmp = column_std(gene3c.X)
    tmp = rankdata(-tmp)
    gene3c.var['highly_variable'] = (tmp < 2000)
    gene3c = gene3c[:, gene3c.var['highly_variable']].copy()

    # Create sparse matrix for output
    gene3c_mtx = csr_matrix(gene3c.X)
    gene3c_cells = gene3c.obs.index
    gene3c_genes = gene3c.var.index

    # Save barcodes and genes
    pd.DataFrame(gene3c_cells).to_csv(f'hicluster/imputed_matrix/{args.res}kb_resolution/genescore/geneimputescore.barcodes.tsv', sep='\t')
    pd.DataFrame(gene3c_genes).to_csv(f'hicluster/imputed_matrix/{args.res}kb_resolution/genescore/geneimputescore.genes.tsv', sep='\t')

    # Write the sparse matrix to a file
    output_mtx_path = f'hicluster/imputed_matrix/{args.res}kb_resolution/genescore/geneimputescore.mtx'
    mmwrite(output_mtx_path, gene3c_mtx, precision=2)

if __name__ == "__main__":
    main()
