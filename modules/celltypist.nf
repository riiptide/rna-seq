process CELLTYPIST {
    tag "${meta.id}"
    label 'analysis'
    publishDir { "${params.outdir}/analysis/${meta.id}" }, mode: 'copy', pattern: '*.{h5ad,csv,json}'
    publishDir "${params.outdir}/analysis", mode: 'copy', pattern: 'annotated.h5ad', enabled: { meta.id == params.sample_id }
    publishDir params.outdir, mode: 'copy', pattern: 'annotated.h5ad', enabled: { meta.id == params.sample_id }
    publishDir { "${params.outdir}/figures/${meta.id}/annotation" }, mode: 'copy', pattern: 'figures/*', saveAs: { name -> name.tokenize('/').last() }
    input:
    tuple val(meta), path(clustered), path(annotation_counts)
    val options
    path code
    output:
    tuple val(meta), path('annotated.h5ad'), emit: annotated
    tuple val(meta), path('annotation_summary.json'), emit: summary
    path '*.csv', emit: tables
    path 'figures/*', emit: figures
    script:
    """
    python celltypist_annotate.py --input '${clustered}' --annotation-counts '${annotation_counts}' ${options}
    """
}
