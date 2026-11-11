process MULTIQC {
    tag 'all_samples'
    label 'light'
    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true
    input:
    path reports, stageAs: 'reports??/*'
    output:
    path 'multiqc_report.html', emit: html
    path 'multiqc_report_data', emit: data
    script:
    """
    multiqc . --force --filename multiqc_report.html
    """
}
