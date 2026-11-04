process SCANPY_CLUSTER {
    tag "${meta.id}"
    label 'analysis'
    publishDir { "${params.outdir}/analysis/${meta.id}" }, mode: 'copy', pattern: '*.{h5ad,csv,json}'
    publishDir "${params.outdir}/analysis", mode: 'copy', pattern: 'markers.csv', enabled: { meta.id == params.sample_id }
    publishDir { "${params.outdir}/figures/${meta.id}/cluster" }, mode: 'copy', pattern: 'figures/*', saveAs: { name -> name.tokenize('/').last() }
    input:
    tuple val(meta), path(filtered)
    val options
    path code
    output:
    tuple val(meta), path('clustered.h5ad'), emit: clustered
    tuple val(meta), path('markers.csv'), emit: markers
    path 'markers_top.csv', emit: top_markers
    path 'pca_variance.csv', emit: variance
    path 'cluster_summary.json', emit: summary
    path 'figures/*', emit: figures
    script:
    """
    python scanpy_cluster.py --input '${filtered}' ${options}
    """
}
