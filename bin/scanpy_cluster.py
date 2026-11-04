#!/usr/bin/env python3
"""Normalize, select HVGs, reduce dimensions, cluster, and rank marker genes."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from pipeline_utils import record, save_h5ad, save_json, seed_all


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True)
    p.add_argument('--output', default='clustered.h5ad')
    for name, default in [('normalization-target', 10000), ('n-top-genes', 2000),
                          ('pca-components', 50), ('n-pcs', 20), ('n-neighbors', 15),
                          ('leiden-iterations', 2), ('marker-top-n', 4), ('seed', 42),
                          ('plot-dpi', 160)]:
        p.add_argument('--'+name, type=int, default=default)
    p.add_argument('--hvg-flavor', default='seurat')
    p.add_argument('--pca-scale-clip', type=float, default=10)
    p.add_argument('--leiden-resolution', type=float, default=0.6)
    p.add_argument('--umap-min-dist', type=float, default=0.3)
    p.add_argument('--marker-method', default='wilcoxon')
    a = p.parse_args()
    seed_all(a.seed)
    Path('figures').mkdir(exist_ok=True)
    adata = sc.read_h5ad(a.input)
    sc.pp.normalize_total(adata, target_sum=a.normalization_target)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    sc.pp.highly_variable_genes(adata, flavor=a.hvg_flavor, n_top_genes=min(a.n_top_genes, adata.n_vars))
    adata.var['highly_variable'] &= ~adata.var.mt
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig('figures/hvgs.png', dpi=a.plot_dpi, bbox_inches='tight')
    plt.close('all')
    hvg = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(hvg, max_value=a.pca_scale_clip)
    components = min(a.pca_components, hvg.n_obs-1, hvg.n_vars-1)
    if components < 2:
        raise ValueError('Not enough cells/HVGs for PCA')
    sc.tl.pca(hvg, n_comps=components, svd_solver='arpack', random_state=a.seed)
    adata.obsm['X_pca'] = hvg.obsm['X_pca'].copy()
    adata.uns['pca'] = hvg.uns['pca'].copy()
    adata.varm['PCs'] = np.zeros((adata.n_vars, components), dtype=np.float32)
    adata.varm['PCs'][adata.var.highly_variable.to_numpy(), :] = hvg.varm['PCs']
    pcs = min(a.n_pcs, components)
    neighbors = min(a.n_neighbors, adata.n_obs-1)
    variance = pd.DataFrame({'PC': np.arange(1, components+1),
                            'variance_ratio': adata.uns['pca']['variance_ratio']})
    variance.to_csv('pca_variance.csv', index=False)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
    for ax in axes:
        ax.plot(variance.PC, variance.variance_ratio, 'o-', markersize=3)
        ax.axvline(pcs, color='crimson', linestyle='--')
        ax.set(xlabel='Principal component', ylabel='Explained variance')
    axes[1].set_yscale('log')
    fig.savefig('figures/pca_elbow.png', dpi=a.plot_dpi)
    plt.close(fig)
    sc.pp.neighbors(adata, n_neighbors=neighbors, n_pcs=pcs, use_rep='X_pca', random_state=a.seed)
    sc.tl.leiden(adata, resolution=a.leiden_resolution, random_state=a.seed, flavor='igraph',
                 directed=False, n_iterations=a.leiden_iterations, key_added='leiden')
    sc.tl.umap(adata, min_dist=a.umap_min_dist, random_state=a.seed)
    if adata.obs.leiden.nunique() < 2 or adata.obs.leiden.value_counts().min() < 2:
        raise ValueError('Marker testing requires at least two clusters with at least two cells each')
    sc.tl.rank_genes_groups(adata, groupby='leiden', method=a.marker_method, use_raw=True,
                           tie_correct=True, corr_method='benjamini-hochberg', pts=True)
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv('markers.csv', index=False)
    markers[markers.scores > 0].groupby('group', observed=True).head(a.marker_top_n).to_csv('markers_top.csv', index=False)
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=a.marker_top_n, standard_scale='var', dendrogram=False, show=False)
    plt.savefig('figures/cluster_markers.png', dpi=a.plot_dpi, bbox_inches='tight')
    plt.close('all')
    colors = ['leiden', 'total_counts', 'n_genes_by_counts', 'pct_counts_mt']
    if 'doublet_score' in adata.obs:
        colors += ['doublet_score']
    sc.pl.umap(adata, color=colors, ncols=2, show=False)
    plt.savefig('figures/umap_clusters_qc.png', dpi=a.plot_dpi, bbox_inches='tight')
    plt.close('all')
    record(adata, 'cluster', a)
    adata.uns['effective_dimensions'] = {'pca_components': components, 'neighbor_pcs': pcs,
                                       'neighbors': neighbors, 'hvgs': hvg.n_vars}
    save_h5ad(adata, a.output)
    save_json('cluster_summary.json', {'cells': adata.n_obs, 'clusters': int(adata.obs.leiden.nunique())})


if __name__ == '__main__':
    main()
