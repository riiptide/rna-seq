#!/usr/bin/env python3
"""Annotate with a checksum-pinned, offline CellTypist model and graph majority voting."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import celltypist
from celltypist.models import Model
from pipeline_utils import checksum, record, save_h5ad, save_json, seed_all


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True)
    p.add_argument('--annotation-counts', required=True)
    p.add_argument('--output', default='annotated.h5ad')
    p.add_argument('--model', default='/opt/conda/models/Immune_All_Low.pkl')
    p.add_argument('--model-sha256', default='290874d35dac039d4c9218c343fde4aac1077709b72a331ce7266f6828c36502')
    p.add_argument('--normalization-target', type=int, default=10000)
    p.add_argument('--overcluster-resolution', type=float, default=5.0)
    p.add_argument('--min-prop', type=float, default=0.5)
    p.add_argument('--leiden-iterations', type=int, default=2)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--plot-dpi', type=int, default=160)
    p.add_argument('--pbmc-markers', default='CD3E,IL7R,CD8A,MS4A1,CD79A,GNLY,NKG7,CD14,LYZ,FCGR3A,FCER1A,CST3')
    p.add_argument('--review-markers', default='CD3D,CD3E,TRAC,TRDC,TRGC2,CD8B,KLRD1,FCER1G,TYROBP,SLC4A10,KLRB1,LILRA4,CLEC4C,TCF4,GZMB,CD1C,CCR7,TCF7,LEF1,S100A8,S100A9')
    a = p.parse_args()
    if a.normalization_target != 10000:
        p.error('CellTypist requires log1p of 10,000 counts per cell')
    seed_all(a.seed)
    Path('figures').mkdir(exist_ok=True)
    if checksum(a.model) != a.model_sha256:
        raise ValueError('CellTypist model checksum mismatch')
    adata = sc.read_h5ad(a.input)
    annotation = sc.read_h5ad(a.annotation_counts)
    if not adata.obs_names.equals(annotation.obs_names):
        raise ValueError('Cluster and classifier barcode order differs')
    sc.pp.normalize_total(annotation, target_sum=a.normalization_target)
    sc.pp.log1p(annotation)
    model = Model.load(a.model)
    overlap = np.intersect1d(annotation.var_names, model.features)
    if not len(overlap):
        raise ValueError('No model features overlap; check gene symbols')
    sc.tl.leiden(adata, resolution=a.overcluster_resolution, random_state=a.seed,
                 flavor='igraph', directed=False, n_iterations=a.leiden_iterations,
                 key_added='celltypist_overclustering')
    prediction = celltypist.annotate(annotation, model=model, majority_voting=True,
        over_clustering=adata.obs.celltypist_overclustering, min_prop=a.min_prop)
    labels = prediction.predicted_labels.reindex(adata.obs_names)
    adata.obs['celltypist_predicted'] = labels.predicted_labels.astype('category')
    adata.obs['cell_type'] = labels.majority_voting.astype('category')
    adata.obs['celltypist_best_score'] = prediction.probability_matrix.reindex(adata.obs_names).max(axis=1)
    adata.obs['celltypist_vote_changed'] = adata.obs.celltypist_predicted.astype(str) != adata.obs.cell_type.astype(str)
    labels.to_csv('celltypist_labels.csv')
    pd.crosstab(adata.obs.leiden, adata.obs.cell_type).to_csv('cluster_celltypist_crosstab.csv')
    pd.crosstab(adata.obs.leiden, adata.obs.celltypist_predicted).to_csv('cluster_raw_celltypist_crosstab.csv')
    rows = []
    for cluster in adata.obs.leiden.cat.categories:
        cells = adata.obs[adata.obs.leiden == cluster]
        counts = cells.cell_type.value_counts()
        rows.append({'cluster': cluster, 'cells': len(cells), 'dominant_celltypist_label': str(counts.index[0]),
                     'dominant_label_fraction': float(counts.iloc[0]/len(cells)),
                     'doublets_flagged': int(cells.predicted_doublet.sum()) if 'predicted_doublet' in cells else 0})
    pd.DataFrame(rows).to_csv('cluster_summary.csv', index=False)
    sc.pl.umap(adata, color=['leiden', 'cell_type'], ncols=1, legend_fontsize=8, show=False)
    plt.savefig('figures/umap_celltypes.png', dpi=a.plot_dpi, bbox_inches='tight')
    plt.close('all')
    for name, requested in [('pbmc', a.pbmc_markers), ('followup', a.review_markers)]:
        present = [g for g in requested.split(',') if g in adata.raw.var_names]
        if not present:
            continue
        sc.pl.dotplot(adata, present, groupby='leiden', use_raw=True, standard_scale='var', show=False)
        plt.savefig(f'figures/{name}_markers.png', dpi=a.plot_dpi, bbox_inches='tight')
        plt.close('all')
        values = adata.raw[:, present].X.toarray()
        pd.DataFrame(values, index=adata.obs_names, columns=present).groupby(adata.obs.leiden, observed=True).mean().to_csv(f'{name}_marker_means.csv')
        pd.DataFrame(values > 0, index=adata.obs_names, columns=present).groupby(adata.obs.leiden, observed=True).mean().to_csv(f'{name}_marker_fractions.csv')
    record(adata, 'celltypist', a)
    adata.uns['celltypist_model_sha256'] = a.model_sha256
    if adata.obs[['leiden', 'cell_type']].isna().any().any() or not np.isfinite(adata.obsm['X_umap']).all():
        raise ValueError('Missing labels or nonfinite UMAP')
    save_h5ad(adata, a.output)
    save_json('annotation_summary.json', {'cells': adata.n_obs, 'clusters': int(adata.obs.leiden.nunique()),
        'labels': int(adata.obs.cell_type.nunique()), 'model_feature_overlap': len(overlap),
        'model_features': len(model.features), 'output_sha256': checksum(a.output)})


if __name__ == '__main__':
    main()
