#!/usr/bin/env python3
"""Calculate UMI QC, filter barcodes/genes, and optionally flag/remove Scrublet doublets."""
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
    p.add_argument('--sample-id', default='sample')
    p.add_argument('--output', default='filtered.h5ad')
    p.add_argument('--annotation-output', default='annotation_counts.h5ad')
    for name, default in [('min-genes', 200), ('max-genes', 5000), ('min-counts', 1000),
                          ('max-counts', 30000), ('gene-min-cells', 3), ('seed', 42),
                          ('doublet-n-pcs', 30), ('plot-bins', 60), ('plot-dpi', 160)]:
        p.add_argument('--'+name, type=int, default=default)
    p.add_argument('--max-pct-mt', type=float, default=20)
    p.add_argument('--mt-prefix', default='MT-')
    p.add_argument('--run-doublets', action=argparse.BooleanOptionalAction, default=True)
    p.add_argument('--remove-doublets', action=argparse.BooleanOptionalAction, default=False)
    p.add_argument('--doublet-expected-rate', type=float, default=0.01)
    p.add_argument('--doublet-sim-ratio', type=float, default=2.0)
    p.add_argument('--doublet-threshold', type=float, default=0.0, help='0 chooses Scrublet automatic threshold')
    a = p.parse_args()
    if a.remove_doublets and not a.run_doublets:
        p.error('--remove-doublets requires --run-doublets')
    seed_all(a.seed)
    Path('figures').mkdir(exist_ok=True)
    adata = sc.read_h5ad(a.input)
    adata.var['mt'] = adata.var.gene_symbol.str.startswith(a.mt_prefix)
    if not adata.var.mt.any():
        raise ValueError('No mitochondrial genes mapped: inspect t2g/symbols')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    columns = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    qc = adata.obs[columns].copy()
    criteria = pd.DataFrame({
        'genes_in_range': qc.n_genes_by_counts.between(a.min_genes, a.max_genes),
        'counts_in_range': qc.total_counts.between(a.min_counts, a.max_counts),
        'mitochondrial_in_range': qc.pct_counts_mt.le(a.max_pct_mt)}, index=qc.index)
    qc['pass_qc'] = criteria.all(axis=1)
    qc.join(criteria).to_csv('barcode_qc.csv.gz')
    candidates = qc[qc.n_genes_by_counts >= a.min_genes]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout='constrained')
    for ax, col, limits in zip(axes, columns, [[a.min_genes, a.max_genes],
                              [a.min_counts, a.max_counts], [a.max_pct_mt]]):
        ax.hist(candidates[col], bins=a.plot_bins)
        for limit in limits:
            ax.axvline(limit, color='crimson', linestyle='--')
        ax.set(xlabel=col, ylabel='Candidate barcodes')
    fig.savefig('figures/qc_before.png', dpi=a.plot_dpi)
    plt.close(fig)
    ranked = np.sort(qc.total_counts.to_numpy())[::-1]
    ranked = ranked[ranked > 0]
    fig, ax = plt.subplots()
    ax.loglog(np.arange(1, len(ranked)+1), ranked)
    ax.axhline(a.min_counts, color='crimson', linestyle='--')
    ax.set(xlabel='Barcode rank', ylabel='UMIs')
    fig.savefig('figures/barcode_rank.png', dpi=a.plot_dpi)
    plt.close(fig)
    adata = adata[qc.pass_qc].copy()
    if adata.n_obs < 4:
        raise ValueError(f'Only {adata.n_obs} cells passed QC; inspect thresholds/read depth')
    summary = {'input_barcodes': len(qc), 'cells_after_thresholds': adata.n_obs,
               'doublets_run': a.run_doublets, 'doublets_flagged': 0, 'doublets_removed': 0}
    if a.run_doublets:
        # Scrublet takes raw counts and handles its own feature selection/normalization.
        sc.pp.scrublet(adata, expected_doublet_rate=a.doublet_expected_rate,
                      sim_doublet_ratio=a.doublet_sim_ratio, random_state=a.seed,
                      n_prin_comps=min(a.doublet_n_pcs, adata.n_obs-1),
                      threshold=a.doublet_threshold or None, use_approx_neighbors=False)
        if adata.obs.predicted_doublet.isna().any():
            raise ValueError('Scrublet did not produce valid calls')
        summary['doublets_flagged'] = int(adata.obs.predicted_doublet.sum())
        summary['doublet_threshold'] = float(adata.uns['scrublet']['threshold'])
        adata.obs[['doublet_score', 'predicted_doublet']].to_csv('doublets.csv')
        fig, ax = plt.subplots()
        ax.hist(adata.obs.doublet_score, bins=a.plot_bins)
        ax.axvline(summary['doublet_threshold'], color='crimson', linestyle='--')
        ax.set(xlabel='Scrublet score', ylabel='Cells')
        fig.savefig('figures/doublet_scores.png', dpi=a.plot_dpi)
        plt.close(fig)
        if a.remove_doublets:
            summary['doublets_removed'] = summary['doublets_flagged']
            adata = adata[~adata.obs.predicted_doublet].copy()
    record(adata, 'qc', a)
    adata.uns['qc_summary'] = summary
    # Preserve all mapped features for the classifier, separately from the gene filter.
    save_h5ad(adata, a.annotation_output)
    sc.pp.filter_genes(adata, min_cells=a.gene_min_cells)
    adata.layers['counts'] = adata.X.copy()
    summary.update(cells_after_qc=adata.n_obs, genes_after_qc=adata.n_vars)
    adata.uns['qc_summary'] = summary
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout='constrained')
    for ax, col in zip(axes, columns):
        ax.hist(adata.obs[col], bins=a.plot_bins)
        ax.set(xlabel=col, ylabel='Retained cells')
    fig.savefig('figures/qc_after.png', dpi=a.plot_dpi)
    plt.close(fig)
    save_h5ad(adata, a.output)
    save_json('qc_summary.json', summary)
    custom = {'id': 'cell_qc', 'section_name': 'Cell QC and doublets', 'plot_type': 'table',
              'data': {a.sample_id: {key: value for key, value in summary.items() if not isinstance(value, bool)}}}
    save_json('cell_qc_mqc.json', custom)
    print(summary)


if __name__ == '__main__':
    main()
