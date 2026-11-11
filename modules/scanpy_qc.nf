process SCANPY_QC {
    tag "${meta.id}"
    label 'analysis'
    publishDir { "${params.outdir}/figures/${meta.id}/qc" }, mode: 'copy', overwrite: true, pattern: 'figures/*', saveAs: { name -> name.tokenize('/').last() }
    publishDir { "${params.outdir}/analysis/${meta.id}/qc" }, mode: 'copy', overwrite: true, pattern: '*.{json,csv,csv.gz}'
    input:
    tuple val(meta), path(counts)
    val options
    path code
    output:
    tuple val(meta), path('filtered.h5ad'), emit: filtered
    tuple val(meta), path('annotation_counts.h5ad'), emit: annotation
    tuple val(meta), path('cell_qc_mqc.json'), emit: multiqc
    path 'qc_summary.json', emit: summary
    path 'barcode_qc.csv.gz', emit: barcode_qc
    path 'doublets.csv', optional: true, emit: doublets
    path 'figures/*', emit: figures
    script:
    """
    python scanpy_qc.py --input '${counts}' --sample-id '${meta.id}' ${options}
    """
}
