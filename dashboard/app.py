"""Interactive exploration of the saved Nextflow PBMC analysis."""
import argparse
import sys
from html import escape
from pathlib import Path
import streamlit as st
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dashboard.charts import colors_for, qc_histogram, umap_figure
from dashboard.data import cell_frame, cluster_order, gene_expression, qc_thresholds, read_analysis, read_markers, summary

CHART_CONFIG = {'displaylogo': False, 'responsive': True, 'scrollZoom': False,
                'modeBarButtonsToRemove': ['select2d', 'lasso2d']}


@st.cache_resource(show_spinner='Loading the annotated cells…', max_entries=2)
def load_analysis(path, fingerprint):
    return read_analysis(path)


@st.cache_data(show_spinner=False, max_entries=2)
def load_markers(path, fingerprint):
    return read_markers(path)


def fingerprint(path):
    stat = path.stat()
    return (stat.st_mtime_ns, stat.st_size)


def legend(frame, column, colors):
    counts = frame[column].value_counts()
    parts = []
    for label, color in colors.items():
        name = f'Cluster {label}' if column == 'leiden' else label
        parts.append(f'<span class="pbmc-swatch" style="background:{color}"></span>'
                     f'<span>{escape(name)}</span><span class="pbmc-count">{counts[label]:,}</span>')
    st.markdown('<div class="pbmc-legend">'+''.join(parts)+'</div>', unsafe_allow_html=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=ROOT/'results/analysis/annotated.h5ad')
    parser.add_argument('--markers', type=Path, default=ROOT/'results/analysis/markers.csv')
    args = parser.parse_args()
    st.set_page_config(page_title='PBMC explorer', page_icon='🧬', layout='wide')
    st.markdown('<style>'+Path(__file__).with_name('style.css').read_text()+'</style>', unsafe_allow_html=True)
    st.markdown('<p class="pbmc-kicker">Single-cell RNA sequencing · 10x Genomics PBMC</p>', unsafe_allow_html=True)
    st.title('PBMC explorer')
    st.caption('Explore immune-cell identities, gene expression, and the quality checks behind this sample.')
    try:
        data_path = args.data.resolve()
        adata = load_analysis(str(data_path), fingerprint(data_path))
    except (OSError, ValueError, KeyError) as error:
        st.error(f'Could not load the annotated analysis. {error}')
        st.info('Generate the full pipeline output from the repository root, or pass --data with an annotated .h5ad path.')
        st.code('nextflow run main.nf -profile docker', language='bash')
        st.stop()
    frame = cell_frame(adata)
    totals = summary(frame)
    metric_columns = st.columns(4)
    for column, label, value in zip(metric_columns,
            ['Cells passing QC', 'Leiden clusters', 'Median genes / cell', 'Predicted cell types'],
            [totals['cells'], totals['clusters'], totals['median_genes'], totals['types']]):
        column.metric(label, f'{value:,.0f}')
    if totals['unassigned']:
        st.caption(f"{totals['unassigned']:,} cells have heterogeneous or unassigned labels and are included in the cell total.")
    umap_tab, qc_tab, marker_tab = st.tabs(['UMAP & expression', 'Quality control', 'Marker genes'])
    with umap_tab:
        with st.container(border=True):
            mode = st.radio('Color cells by', ['Cell type', 'Leiden cluster', 'Gene expression'], horizontal=True)
            if mode == 'Gene expression':
                genes = sorted(adata.raw.var_names.astype(str))
                default_gene = genes.index('MS4A1') if 'MS4A1' in genes else 0
                gene = st.selectbox('Search gene', genes, index=default_gene, help='Type a symbol to search the genes saved in the analysis.')
                values = gene_expression(adata, gene)
                plotted = frame.assign(expression=values)
                st.plotly_chart(umap_figure(plotted, 'expression', gene=gene), width='stretch', config=CHART_CONFIG, key='gene_umap')
                detected = int((values > 0).sum())
                st.caption(f'{gene} is detected in {detected:,} of {len(frame):,} cells ({detected/len(frame):.1%}). Color shows log1p expression after normalization to 10,000 UMIs per cell.')
            else:
                column = 'cell_type' if mode == 'Cell type' else 'leiden'
                colors = colors_for(frame, column)
                plot_column, legend_column = st.columns([3.2, 1.2], gap='large')
                with plot_column:
                    st.plotly_chart(umap_figure(frame, column, colors=colors), width='stretch', config=CHART_CONFIG, key='category_umap')
                with legend_column:
                    st.markdown('**Cell types**' if column == 'cell_type' else '**Leiden clusters**')
                    legend(frame, column, colors)
                st.caption('Each dot is one cell. Hover for its barcode, labels, and QC metrics. Drag to pan and use the chart toolbar to zoom or reset.')
    with qc_tab:
        st.subheader('Quality of retained cells')
        st.caption('Histograms show the cells saved after pipeline filtering. Dashed lines mark the recorded cutoffs.')
        limits = qc_thresholds(adata)
        for column, metric, label in zip(st.columns(3),
                ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                ['Detected genes per cell', 'UMIs per cell', 'Mitochondrial UMIs (%)']):
            with column:
                st.plotly_chart(qc_histogram(frame[metric], label, limits[metric]), width='stretch', config=CHART_CONFIG, key=metric)
                if any(value is None for value in limits[metric]):
                    st.caption('Some filter thresholds were absent from this file’s metadata.')
        if 'predicted_doublet' in adata.obs:
            flagged = int(adata.obs.predicted_doublet.sum())
            st.info(f'{flagged:,} displayed cells carry a Scrublet doublet flag. The saved run records whether flagged cells were removed in its QC settings.')
        with st.expander('Recorded QC settings'):
            st.json(dict(adata.uns.get('qc_params', adata.uns.get('qc_thresholds', {}))))
    with marker_tab:
        st.subheader('Genes enriched in each cluster')
        st.caption('Positive Wilcoxon scores rank genes expressed more strongly in a cluster. Adjusted p-values use the Benjamini–Hochberg correction.')
        try:
            marker_path = args.markers.resolve()
            markers = load_markers(str(marker_path), fingerprint(marker_path))
        except (OSError, ValueError) as error:
            st.warning(f'Could not load the marker table. {error}')
        else:
            cluster = st.selectbox('Leiden cluster', cluster_order(frame.leiden))
            cluster_cells = frame.loc[frame.leiden == cluster]
            leader = cluster_cells.cell_type.value_counts().index[0]
            proportion = (cluster_cells.cell_type == leader).mean()
            st.caption(f'Cluster {cluster} · {len(cluster_cells):,} cells · leading label {leader} ({proportion:.0%})')
            selected = markers.loc[(markers.group == cluster) & (markers.scores > 0)].sort_values('scores', ascending=False, kind='stable')
            limit = st.selectbox('Rows to display', [10, 25, 50, 100], index=1)
            columns = ['names', 'scores', 'logfoldchanges', 'pvals_adj']
            if 'pct_nz_group' in selected.columns:
                columns.append('pct_nz_group')
            shown = selected[columns].head(limit).rename(columns={
                'names': 'Gene', 'scores': 'Wilcoxon score', 'logfoldchanges': 'Log2 fold change',
                'pvals_adj': 'Adjusted p-value', 'pct_nz_group': 'Fraction detected'})
            if shown.empty:
                st.info('This cluster has no positive marker scores in the supplied table.')
            st.dataframe(shown, width='stretch', hide_index=True, column_config={
                'Wilcoxon score': st.column_config.NumberColumn(format='%.2f'),
                'Log2 fold change': st.column_config.NumberColumn(format='%.2f'),
                'Adjusted p-value': st.column_config.NumberColumn(format='%.2e'),
                'Fraction detected': st.column_config.NumberColumn(format='%.2f')})
            st.download_button('Download cluster markers', selected.to_csv(index=False),
                               file_name=f'cluster_{cluster}_markers.csv', mime='text/csv')
    with st.expander('Sample and analysis details'):
        sample = adata.uns.get('qc_params', {}).get('sample_id', data_path.stem)
        st.write(f'Sample: {sample}')
        st.write(f'Annotated matrix: {data_path.name} · {adata.n_obs:,} cells × {adata.n_vars:,} genes')
        st.write('Cell types are CellTypist Immune_All_Low predictions with majority voting. Cluster identities were reviewed against PBMC marker genes during analysis.')
        st.json(dict(adata.uns.get('software_versions', {})))


if __name__ == '__main__':
    main()
