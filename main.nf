include { FASTQC } from './modules/fastqc'
include { KB_REF } from './modules/kb_ref'
include { KB_COUNT } from './modules/kb_count'
include { TO_ANNDATA } from './modules/to_anndata'
include { SCANPY_QC } from './modules/scanpy_qc'
include { SCANPY_CLUSTER } from './modules/scanpy_cluster'
include { CELLTYPIST } from './modules/celltypist'
include { MULTIQC } from './modules/multiqc'

// Quote CLI values; biological parameters are supplied explicitly to each module.
def cli(Map options) {
    return options.collect { name, value -> "--${name} '${value.toString().replace("'", "'\"'\"'")}'" }.join(' ')
}

workflow {
    def index_path = file(params.kb_index ?: params.index)
    def t2g_path = file(params.t2g)
    def reference = channel.value(tuple(index_path, t2g_path))
    if (!index_path.exists() || !t2g_path.exists()) {
        if (params.kb_index) {
            error "Explicit --kb_index/--t2g pair does not exist"
        }
        KB_REF(file("${projectDir}/docs/downloads.sha256"), file("${projectDir}/docs/inputs.sha256"),
               params.reference_url, params.download_retries, file("${projectDir}/bin/reference.py"))
        reference = KB_REF.out.reference.first()
    }
    // groupTuple keeps R1/R2 paired in lane order, and counts all lanes together.
    reads = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2) {
                error 'Samplesheet needs sample, fastq_1, fastq_2 columns'
            }
            if (!(row.sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
                error "Unsafe sample identifier: ${row.sample}"
            }
            tuple(row.sample, [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)])
        }
        .groupTuple()
        .map { id, lanes -> tuple([id: id], lanes.flatten()) }
    FASTQC(reads, params.fastqc_memory_mb)
    KB_COUNT(reads, reference, file("${projectDir}/docs/inputs.sha256"),
             "-x '${params.chemistry}' --workflow '${params.workflow}' -m '${params.sort_memory}'",
             [file("${projectDir}/bin/reference.py"), file("${projectDir}/bin/inspect_counts.py")])
    TO_ANNDATA(KB_COUNT.out.h5ad, reference.map { _idx, t2g -> t2g },
               [file("${projectDir}/bin/to_anndata.py"), file("${projectDir}/bin/pipeline_utils.py")])
    def qc_options = cli([
        'min-genes': params.qc_min_genes, 'max-genes': params.qc_max_genes,
        'min-counts': params.qc_min_counts, 'max-counts': params.qc_max_counts,
        'max-pct-mt': params.qc_max_pct_mt, 'gene-min-cells': params.qc_gene_min_cells,
        'mt-prefix': params.mt_prefix, 'seed': params.seed,
        'doublet-expected-rate': params.doublet_expected_rate, 'doublet-sim-ratio': params.doublet_sim_ratio,
        'doublet-threshold': params.doublet_threshold, 'doublet-n-pcs': params.doublet_n_pcs,
        'plot-bins': params.plot_bins, 'plot-dpi': params.plot_dpi]) +
        (params.run_doublets ? ' --run-doublets' : ' --no-run-doublets') +
        (params.remove_doublets ? ' --remove-doublets' : ' --no-remove-doublets')
    SCANPY_QC(TO_ANNDATA.out.counts, qc_options,
              [file("${projectDir}/bin/scanpy_qc.py"), file("${projectDir}/bin/pipeline_utils.py")])
    SCANPY_CLUSTER(SCANPY_QC.out.filtered, cli([
        'normalization-target': params.normalization_target, 'n-top-genes': params.hvg_n_top_genes,
        'hvg-flavor': params.hvg_flavor, 'pca-components': params.pca_components,
        'pca-scale-clip': params.pca_scale_clip, 'n-pcs': params.neighbors_n_pcs,
        'n-neighbors': params.neighbors_n_neighbors, 'leiden-resolution': params.leiden_resolution,
        'leiden-iterations': params.leiden_iterations, 'umap-min-dist': params.umap_min_dist,
        'marker-method': params.marker_method, 'marker-top-n': params.marker_top_n,
        'seed': params.seed, 'plot-dpi': params.plot_dpi]),
        [file("${projectDir}/bin/scanpy_cluster.py"), file("${projectDir}/bin/pipeline_utils.py")])
    CELLTYPIST(SCANPY_CLUSTER.out.clustered.join(SCANPY_QC.out.annotation), cli([
        'model': params.celltypist_model, 'model-sha256': params.celltypist_model_sha256,
        'normalization-target': params.normalization_target,
        'overcluster-resolution': params.celltypist_overcluster_resolution,
        'min-prop': params.celltypist_min_prop, 'leiden-iterations': params.leiden_iterations,
        'seed': params.seed, 'plot-dpi': params.plot_dpi,
        'pbmc-markers': params.pbmc_markers, 'review-markers': params.review_markers]),
        [file("${projectDir}/bin/celltypist_annotate.py"), file("${projectDir}/bin/pipeline_utils.py")])
    // collect waits until all samples finish, so MultiQC is the final aggregation.
    reports = FASTQC.out.zip.map { _meta, files -> files }
        .mix(SCANPY_QC.out.multiqc.map { _meta, report -> report })
        .mix(CELLTYPIST.out.summary.map { _meta, report -> report })
        .collect()
    MULTIQC(reports)
}
