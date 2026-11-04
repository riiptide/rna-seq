process KB_COUNT {
    tag "${meta.id}"
    label 'count'
    input:
    tuple val(meta), path(reads)
    tuple path(index), path(t2g)
    path checksums
    val options
    path code
    output:
    tuple val(meta), path('kb/counts_unfiltered/adata.h5ad'), emit: h5ad
    tuple val(meta), path('kb/run_info.json'), emit: info
    path 'count_validation.json', emit: validation
    script:
    def files = reads.collect { read -> "'${read}'" }.join(' ')
    """
    python reference.py --index '${index}' --t2g '${t2g}' --inputs-manifest '${checksums}'
    kb count -i '${index}' -g '${t2g}' -t ${task.cpus} ${options} --h5ad -o kb ${files}
    python inspect_counts.py kb/counts_unfiltered/adata.h5ad > count_validation.json
    """
}
