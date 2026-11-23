"""Plotly figures with stable colors, cell hover details, and recorded QC boundaries."""
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from dashboard.data import cluster_order

PALETTE = ['#2563eb', '#0f766e', '#9333ea', '#c2410c', '#be185d', '#4d7c0f',
           '#0891b2', '#7c3aed', '#a16207', '#b91c1c', '#047857', '#475569',
           '#db2777', '#4338ca', '#15803d', '#0369a1', '#9a3412', '#6b21a8',
           '#4f46e5', '#a21caf', '#0e7490', '#92400e', '#166534', '#9f1239']


def colors_for(frame, column):
    labels = cluster_order(frame[column]) if column == 'leiden' else sorted(frame[column].unique())
    return {label: ('#64748b' if label == 'Heterogeneous' else PALETTE[i % len(PALETTE)])
            for i, label in enumerate(labels)}


def umap_figure(frame, column, colors=None, gene=None):
    hover = {'barcode': True, 'leiden': True, 'cell_type': True,
             'n_genes_by_counts': ':,', 'total_counts': ':,',
             'pct_counts_mt': ':.1f', 'UMAP 1': False, 'UMAP 2': False}
    names = {'barcode': 'Cell', 'leiden': 'Cluster', 'cell_type': 'Cell type',
             'n_genes_by_counts': 'Detected genes', 'total_counts': 'UMIs',
             'pct_counts_mt': 'Mitochondrial UMIs (%)', 'expression': gene or 'Expression'}
    if gene:
        frame = frame.sort_values('expression', kind='stable')
        hover['expression'] = ':.2f'
    fig = px.scatter(frame, x='UMAP 1', y='UMAP 2', color=column,
                     color_discrete_map=colors, labels=names, hover_data=hover,
                     color_continuous_scale=['#e2e8f0', '#b9d9d2', '#65af9e', '#16816c', '#064e3b'],
                     range_color=(0, max(float(frame.expression.max()), 1e-6)) if gene else None,
                     render_mode='svg')
    fig.update_traces(marker={'size': 6, 'opacity': 0.85, 'line': {'width': 0}})
    fig.update_layout(template='plotly_white', height=510, showlegend=False,
                      margin={'l': 20, 'r': 15, 't': 10, 'b': 35},
                      paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='#f8fafc',
                      font={'family': 'Inter, Arial, sans-serif', 'size': 14, 'color': '#334155'},
                      uirevision='pbmc-umap', dragmode='pan')
    fig.update_xaxes(showgrid=False, zeroline=False, showticklabels=False)
    fig.update_yaxes(showgrid=False, zeroline=False, showticklabels=False, scaleanchor='x', scaleratio=1)
    if gene:
        fig.update_layout(coloraxis_colorbar={'title': 'Log1p<br>expression', 'thickness': 14, 'len': 0.75})
    return fig


def qc_histogram(values, label, limits, bins=40):
    values = np.asarray(values, dtype=float)
    finite_limits = [v for v in limits if v is not None]
    high = max([float(values.max()), *finite_limits, 1]) * 1.08
    frequencies, edges = np.histogram(values, bins=bins, range=(0, high))
    fig = go.Figure(go.Bar(x=(edges[:-1]+edges[1:])/2, y=frequencies,
                          width=np.diff(edges)*0.93, marker_color='#0f766e',
                          hovertemplate='%{x:.1f}<br>%{y} cells<extra></extra>'))
    for index, limit in enumerate(limits):
        if limit is not None:
            name = 'Max' if len(limits) == 1 or index == 1 else 'Min'
            fig.add_vline(x=limit, line_color='#b45309', line_dash='dash', line_width=2,
                          annotation_text=f'{name} {limit:,.0f}',
                          annotation_position='top left' if name == 'Max' else 'top right')
    fig.update_layout(template='plotly_white', height=310, margin={'l': 45, 'r': 15, 't': 38, 'b': 45},
                      xaxis_title=label, yaxis_title='Cells', bargap=0.06,
                      paper_bgcolor='rgba(0,0,0,0)', font={'size': 14, 'color': '#334155'})
    fig.update_xaxes(range=[0, high], showgrid=False)
    fig.update_yaxes(rangemode='tozero', gridcolor='#edf2f7')
    return fig
