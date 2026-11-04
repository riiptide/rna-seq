process TO_ANNDATA {
    tag "${meta.id}"
    label 'analysis'
    publishDir { "${params.outdir}/counts/${meta.id}" }, mode: 'copy'
    input:
    tuple val(meta), path(counts)
    path t2g
    path code
    output:
    tuple val(meta), path('counts.h5ad'), emit: counts
    tuple val(meta), path('mapping_summary.json'), emit: summary
    path 'gene_id_mapping.tsv', emit: mapping
    script:
    """
    python to_anndata.py --input '${counts}' --t2g '${t2g}' --output counts.h5ad
    """
}
