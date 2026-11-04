#!/usr/bin/env python3
"""Map raw kb AnnData gene IDs to symbols, summing duplicate symbols without losing UMIs."""
import argparse
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse as sp
from pipeline_utils import checksum, save_h5ad, save_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='Raw counts_unfiltered/adata.h5ad from kb count')
    parser.add_argument('--t2g', required=True)
    parser.add_argument('--output', default='counts.h5ad')
    args = parser.parse_args()
    source = sc.read_h5ad(args.input)
    assert sp.issparse(source.X) and source.obs_names.is_unique
    assert np.isfinite(source.X.data).all() and (source.X.data >= 0).all()
    assert np.equal(source.X.data, np.floor(source.X.data)).all()
    input_shape = source.shape
    input_total_umis = float(source.X.sum())
    reference = pd.read_csv(args.t2g, sep='\t', header=None, usecols=[0, 1, 2],
                            names=['transcript_id', 'gene_id', 'symbol'], dtype=str,
                            keep_default_na=False)
    assert reference.groupby('gene_id')['symbol'].nunique().max() == 1
    exact = reference.drop_duplicates('gene_id').set_index('gene_id')['symbol']
    base = reference.assign(base_id=reference.gene_id.str.replace(r'\.\d+$', '', regex=True))
    unambiguous = base.groupby('base_id')['symbol'].nunique().eq(1)
    fallback = base[base.base_id.isin(unambiguous[unambiguous].index)].drop_duplicates(
        'base_id').set_index('base_id')['symbol']
    original_ids = pd.Index(source.var_names.astype(str))
    ensembl = original_ids.str.startswith('ENSG')
    symbols = pd.Series(original_ids, index=original_ids, dtype=object)
    symbols.loc[ensembl] = exact.reindex(original_ids[ensembl]).values
    missing = ensembl & symbols.isna().to_numpy()
    symbols.loc[missing] = fallback.reindex(
        original_ids[missing].str.replace(r'\.\d+$', '', regex=True)).values
    symbols = symbols.fillna('').str.strip()
    feature_names = pd.Index(np.where(symbols.eq(''), original_ids, symbols))
    mapping = pd.DataFrame({'original_id': original_ids, 'gene_symbol': symbols.values,
                            'feature_name': feature_names})
    mapping.to_csv('gene_id_mapping.tsv', sep='\t', index=False)
    codes, unique_names = pd.factorize(feature_names, sort=False)
    projection = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)),
                               shape=(source.n_vars, len(unique_names)))
    grouped = mapping.groupby('feature_name', sort=False)
    var = pd.DataFrame(index=pd.Index(unique_names, name='gene'))
    var['original_gene_ids'] = grouped.original_id.agg(';'.join).reindex(var.index).values
    var['gene_symbol'] = grouped.gene_symbol.first().reindex(var.index).values
    var['n_original_ids'] = grouped.size().reindex(var.index).values
    adata = ad.AnnData((source.X @ projection).tocsr().astype(np.float32),
                       obs=source.obs.copy(), var=var)
    assert np.array_equal(np.asarray(source.X.sum(axis=1)).ravel(),
                          np.asarray(adata.X.sum(axis=1)).ravel())
    assert adata.var_names.is_unique
    mapping_summary = {'input_barcodes': input_shape[0], 'original_gene_ids': input_shape[1],
                       'merged_features': adata.n_vars, 'ids_without_symbols': int(symbols.eq('').sum()),
                       'symbol_groups_combining_ids': int((var.n_original_ids > 1).sum())}
    adata.uns['mapping_provenance'] = {'input_sha256': checksum(args.input),
        't2g_sha256': checksum(args.t2g), **mapping_summary}
    save_json('mapping_summary.json', mapping_summary)
    save_h5ad(adata, args.output)


if __name__ == '__main__':
    main()
