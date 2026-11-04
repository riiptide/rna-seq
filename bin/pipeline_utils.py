"""Small shared I/O helpers; analysis parameters belong to each CLI/config."""
import hashlib
import importlib.metadata
import json
import random
from pathlib import Path
import anndata as ad
import numpy as np
from scipy import sparse


def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)


def save_json(path, values):
    Path(path).write_text(json.dumps(values, indent=2, default=str) + '\n')


def checksum(path):
    with open(path, 'rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def save_h5ad(adata, path):
    """Write atomically, then reopen to reject truncated or malformed output."""
    path = Path(path)
    temporary = path.with_suffix('.tmp.h5ad')
    adata.write_h5ad(temporary, compression='gzip')
    check = ad.read_h5ad(temporary)
    if check.shape != adata.shape or not check.obs_names.equals(adata.obs_names):
        raise ValueError('AnnData round-trip changed shape/barcodes')
    if not check.var_names.is_unique or not check.obs_names.is_unique:
        raise ValueError('Duplicate matrix identifiers')
    values = check.X.data if sparse.issparse(check.X) else check.X
    if not np.isfinite(values).all():
        raise ValueError('Nonfinite expression values')
    temporary.replace(path)


def record(adata, step, args):
    adata.uns[step + '_params'] = vars(args).copy()
    packages = ['scanpy', 'anndata', 'celltypist', 'numpy', 'scipy', 'pandas',
                'scikit-learn', 'igraph', 'leidenalg', 'scrublet']
    adata.uns['software_versions'] = {p: importlib.metadata.version(p) for p in packages}
