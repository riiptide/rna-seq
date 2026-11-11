process FASTQC {
    tag "${meta.id}"
    label 'light'
    publishDir { "${params.outdir}/fastqc/${meta.id}" }, mode: 'copy', overwrite: true
    input:
    tuple val(meta), path(reads)
    val memory_mb
    output:
    tuple val(meta), path('*_fastqc.zip'), emit: zip
    tuple val(meta), path('*_fastqc.html'), emit: html
    script:
    def files = reads.collect { read -> "'${read}'" }.join(' ')
    """
    fastqc --threads ${task.cpus} --memory ${memory_mb} ${files}
    """
}
