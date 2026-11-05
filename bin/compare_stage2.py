#!/usr/bin/env python3
"""Compare a full pipeline run to the preserved Stage 2 notebook by barcode overlap."""
import argparse
from pathlib import Path
import json
import anndata as ad
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--notebook', required=True)
    p.add_argument('--pipeline', required=True)
    p.add_argument('--review', default='docs/stage-2-cluster-review.csv')
    p.add_argument('--outdir', required=True)
    a = p.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    notebook = ad.read_h5ad(a.notebook)
    pipeline = ad.read_h5ad(a.pipeline)
    shared = notebook.obs_names.intersection(pipeline.obs_names)
    old = notebook.obs.loc[shared, 'leiden'].astype(str)
    new = pipeline.obs.loc[shared, 'leiden'].astype(str)
    cross = pd.crosstab(new, old)
    cross.to_csv(out/'cluster_overlap.csv')
    ni, oi = linear_sum_assignment(-cross.to_numpy())
    matches = {str(cross.index[i]): str(cross.columns[j]) for i, j in zip(ni, oi)}
    review = pd.read_csv(a.review, dtype={'cluster': str}).set_index('cluster')
    rows = []
    for cluster in pipeline.obs.leiden.cat.categories:
        cells = pipeline.obs[pipeline.obs.leiden == cluster]
        old_cluster = matches.get(str(cluster))
        counts = cells.cell_type.value_counts()
        rows.append({'pipeline_cluster': str(cluster), 'pipeline_cells': len(cells),
            'pipeline_dominant_label': str(counts.index[0]),
            'pipeline_dominant_fraction': float(counts.iloc[0]/len(cells)),
            'notebook_cluster': old_cluster or 'unmatched',
            'notebook_cells': int(review.loc[old_cluster, 'cells']) if old_cluster is not None else 0,
            'notebook_interpretation': review.loc[old_cluster, 'interpretation'] if old_cluster is not None else 'unmatched',
            'shared_cells_in_match': int(cross.loc[str(cluster), old_cluster]) if old_cluster is not None else 0,
            'doublets_flagged': int(cells.predicted_doublet.sum()) if 'predicted_doublet' in cells else 0})
    pd.DataFrame(rows).to_csv(out/'cluster_comparison.csv', index=False)
    report = {'notebook_cells': notebook.n_obs, 'pipeline_cells': pipeline.n_obs,
              'shared_cells': len(shared), 'notebook_only_cells': notebook.n_obs-len(shared),
              'pipeline_only_cells': pipeline.n_obs-len(shared),
              'notebook_clusters': int(notebook.obs.leiden.nunique()),
              'pipeline_clusters': int(pipeline.obs.leiden.nunique()),
              'adjusted_rand_index': adjusted_rand_score(old, new),
              'same_voted_label_fraction': float((notebook.obs.loc[shared, 'cell_type'].astype(str) == pipeline.obs.loc[shared, 'cell_type'].astype(str)).mean()),
              'same_raw_label_fraction': float((notebook.obs.loc[shared, 'celltypist_predicted'].astype(str) == pipeline.obs.loc[shared, 'celltypist_predicted'].astype(str)).mean()),
              'doublets_flagged': int(pipeline.obs.predicted_doublet.sum()) if 'predicted_doublet' in pipeline.obs else 0,
              'doublets_removed': int(pipeline.uns['qc_summary']['doublets_removed'])}
    if notebook.var_names.equals(pipeline.var_names) and notebook.obs_names.equals(pipeline.obs_names):
        report['raw_counts_identical'] = (notebook.layers['counts'] != pipeline.layers['counts']).nnz == 0
    report['cluster_count_within_one'] = abs(report['notebook_clusters']-report['pipeline_clusters']) <= 1
    (out/'comparison.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    if not report['cluster_count_within_one']:
        raise ValueError('Cluster count differs by more than one; inspect graph/platform differences')


if __name__ == '__main__':
    main()
