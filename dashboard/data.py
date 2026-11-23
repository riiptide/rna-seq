"""Read the pipeline's saved analysis and extract the values used by the UI."""
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

QC_COLUMNS = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
UNASSIGNED = {'heterogeneous', 'unassigned', 'unknown'}


def read_analysis(path):
    data = ad.read_h5ad(path)
    required = {'leiden', 'cell_type', *QC_COLUMNS}
    missing = required.difference(data.obs.columns)
    if missing:
        raise ValueError('Missing cell metadata: ' + ', '.join(sorted(missing)))
    if data.n_obs == 0 or 'X_umap' not in data.obsm:
        raise ValueError('The analysis needs cells and a saved UMAP embedding.')
    coordinates = np.asarray(data.obsm['X_umap'])
    if coordinates.shape != (data.n_obs, 2) or not np.isfinite(coordinates).all():
        raise ValueError('The saved UMAP must contain two finite coordinates per cell.')
    if data.obs[list(required)].isna().any().any():
        raise ValueError('The analysis contains missing labels or QC values.')
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError('Cell and gene identifiers must be unique.')
    if data.raw is None or not data.raw.var_names.is_unique:
        raise ValueError('Expected log-normalized expression in .raw, as written by the pipeline.')
    return data


def cell_frame(data):
    frame = data.obs[['leiden', 'cell_type', *QC_COLUMNS]].copy()
    frame['leiden'] = frame.leiden.astype(str)
    frame['cell_type'] = frame.cell_type.astype(str)
    frame['barcode'] = data.obs_names.astype(str)
    frame['UMAP 1'], frame['UMAP 2'] = np.asarray(data.obsm['X_umap']).T
    return frame


def gene_expression(data, gene):
    """Densify just one gene vector, preserving the order of cells in the embedding."""
    if gene not in data.raw.var_names:
        raise ValueError(f'Gene {gene!r} is absent from the saved expression matrix.')
    vector = data.raw[:, gene].X
    values = vector.toarray().ravel() if sparse.issparse(vector) else np.asarray(vector).ravel()
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Expected finite, nonnegative log-normalized gene expression.')
    return values


def qc_thresholds(data):
    values = data.uns.get('qc_params', data.uns.get('qc_thresholds', {}))
    def value(name):
        stored = values.get(name, values.get('qc_' + name))
        return float(stored) if stored is not None else None
    return {
        'n_genes_by_counts': [value('min_genes'), value('max_genes')],
        'total_counts': [value('min_counts'), value('max_counts')],
        'pct_counts_mt': [value('max_pct_mt')],
    }


def read_markers(path):
    frame = pd.read_csv(path, dtype={'group': str, 'names': str})
    required = {'group', 'names', 'scores', 'logfoldchanges', 'pvals_adj'}
    if not required.issubset(frame.columns):
        raise ValueError('The marker table must contain group, names, scores, logfoldchanges and pvals_adj.')
    return frame


def summary(frame):
    labels = set(frame.cell_type)
    unresolved = frame.cell_type.str.lower().isin(UNASSIGNED)
    return {
        'cells': len(frame),
        'clusters': frame.leiden.nunique(),
        'median_genes': float(frame.n_genes_by_counts.median()),
        'types': sum(label.lower() not in UNASSIGNED for label in labels),
        'unassigned': int(unresolved.sum()),
    }


def cluster_order(labels):
    return sorted(set(labels), key=lambda label: (0, int(label)) if str(label).isdigit() else (1, str(label)))
