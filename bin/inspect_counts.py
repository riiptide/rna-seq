#!/usr/bin/env python3
"""Validate a sparse, nonnegative raw UMI matrix and print a small JSON summary."""

import argparse
import json
from pathlib import Path

import anndata
import numpy as np
from scipy import sparse


def inspect(path):
    data = anndata.read_h5ad(path)
    if not sparse.issparse(data.X):
        raise ValueError('Expected sparse counts; refusing to densify the matrix')
    if not data.n_obs or not data.n_vars or not data.X.nnz:
        raise ValueError('Count matrix is empty')
    values = data.X.data
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Counts contain negative or nonfinite values')
    if not (values > 0).any():
        raise ValueError('Count matrix contains no positive counts')
    if not np.equal(values, np.floor(values)).all():
        raise ValueError('Expected integer-valued raw UMI counts')
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError('Barcode and gene identifiers must be unique')
    totals = np.asarray(data.X.sum(axis=1)).ravel()
    return {
        'path': str(path),
        'barcodes_not_called_cells': data.n_obs,
        'reference_genes': data.n_vars,
        'genes_with_counts': int(np.asarray(data.X.sum(axis=0)).ravel().astype(bool).sum()),
        'nonzero_entries': int(data.X.nnz),
        'total_umis': int(totals.sum()),
        'zero_count_barcodes': int((totals == 0).sum()),
        'matrix_format': data.X.format,
        'matrix_dtype': str(data.X.dtype),
        'gene_metadata_columns': list(data.var.columns),
        'example_gene_ids': data.var_names[:3].tolist(),
        'raw_counts_validated': True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('h5ad', type=Path)
    args = parser.parse_args()
    print(json.dumps(inspect(args.h5ad), indent=2))


if __name__ == '__main__':
    main()
