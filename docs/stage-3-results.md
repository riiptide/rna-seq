# Stage 3 validation

The full Docker run processed 66,601,887 read pairs and retained the same 1,107 barcodes and 17,581 genes as the notebook. All 13 Leiden clusters had identical membership (adjusted Rand index 1.0). Every cell had the same raw and majority-voted CellTypist label. Both the test and full unfiltered count matrices matched the Stage 1 matrices element for element.

| Run on Apple Silicon, linux/amd64 containers | Time | Evidence |
|---|---:|---|
| Fresh 1M-pair test with prepared inputs and images | 2m 51s | Seven completed tasks in the initial timed run |
| Full 66.6M-pair run | 17m 40s | Seven completed tasks. FastQC took 17m 31s and ran concurrently with counting, which took 12m 35s |
| Final test with `-resume` | 4.12s | All seven tasks cached. The wall and monotonic timers agreed |

The initial test met the five-minute target before the final output-path fix. Two later fresh runs verified the corrected paths but produced inconsistent clock measurements. One recorded 169s on the monotonic timer and 9m 16s in Nextflow. The last recorded 176.65s monotonic and 3h 21m 21s by UTC, even with idle sleep inhibited. The source of the gaps is unresolved. A consistent sub-five-minute wall time for the final revision remains unverified. Image builds, fixture preparation and reference downloads are separate setup costs. See the retained [timings, traces and logs](stage-3-validation/).

## Cluster comparison

Each row pairs a pipeline cluster with its notebook counterpart by shared cell barcodes. The leading model label summarizes majority-voted predictions within the pipeline cluster. The marker interpretation preserves the notebook's review of mixed states and small groups.

| Pipeline cluster | Notebook cluster | Cells in each | Leading pipeline label | Notebook marker interpretation |
|---:|---:|---:|---|---|
| 0 | 0 | 115 | Classical monocytes | Classical CD14+ monocytes, S100A8/A9-high |
| 1 | 1 | 126 | Naive B cells | Naive B cells |
| 2 | 2 | 76 | Tem/Trm cytotoxic T cells | Cytotoxic T cells — mixed effector/memory-like states |
| 3 | 3 | 174 | Tcm/Naive helper T cells | IL7R+ helper T cells — mixed naive/memory/effector-like states |
| 4 | 4 | 57 | Memory B cells | Memory B cells |
| 5 | 5 | 208 | Classical monocytes | Classical CD14+ monocytes, CPVL/LGALS2-high |
| 6 | 6 | 38 | Tcm/Naive cytotoxic T cells | Naive/central-memory-like CD8 T cells |
| 7 | 7 | 68 | MAIT cells | MAIT-like T cells |
| 8 | 8 | 57 | CD16+ NK cells | NK-like cytotoxic cells — review T-cell-associated RNA |
| 9 | 9 | 32 | Non-classical monocytes | Nonclassical FCGR3A+ monocytes |
| 10 | 10 | 8 | pDC | Plasmacytoid dendritic cells (pDC), small group |
| 11 | 11 | 130 | Tcm/Naive helper T cells | Naive/central-memory-like helper T cells |
| 12 | 12 | 18 | DC2 | Conventional dendritic cells, DC2-like — subtype mixed |

Scrublet flagged 9 of 1,107 cells (0.81%). The default retains those cells with their scores and calls. One flag belonged to C2, four to C5, two to C7 and two to C11. Cluster membership stayed identical under this policy. Excluding the flagged cells would be a separate analysis. The shallow fixture retained 1,109 cells under relaxed QC and yielded four clusters.

## Repository tree

```text
rna-seq/
├── README.md, main.nf, nextflow.config, samplesheet.csv
├── environment.yml, environment-native.yml, environment-stage2.yml
├── notebooks/01_explore.ipynb
├── modules/
│   ├── fastqc.nf, kb_ref.nf, kb_count.nf, to_anndata.nf
│   └── scanpy_qc.nf, scanpy_cluster.nf, celltypist.nf, multiqc.nf
├── containers/
│   ├── build.sh
│   ├── analysis/Dockerfile, analysis/requirements*.txt
│   └── kb/Dockerfile
├── bin/
│   ├── to_anndata.py, scanpy_qc.py, scanpy_cluster.py, celltypist_annotate.py
│   └── pipeline_utils.py, reference.py, inspect_counts.py, validate_fastq_pairs.py, compare_stage2.py
├── tests/make_test_data.sh, tests/samplesheet_test.csv
├── docs/         Walkthroughs, checksum manifests, cluster review and Stage 3 evidence
├── dashboard/, addons/      Reserved for later stages
└── data/, references/, results/, work/      Local data, outputs and cache
```

## Reproduction details

The original Linux kb image carried kallisto 0.52.0 and lacked two runtime libraries. A small derived image supplies the libraries and pins kallisto 0.51.1, matching the Stage 1 count. The reference fallback fetches the checksum-pinned Ensembl 108 archive, and the CellTypist model loads inside its image with networking disabled.

The workflow uses `publishDir` to copy generated files to per-sample directories and create convenience paths for the default sample. Full and test output checks confirmed byte-identical copies of the annotated matrix and marker table at their documented paths. A full cached rerun verified the publishing correction. The Slurm/Apptainer profile awaits a cluster test.

Run commands, parameter rationale and module connections are in the [Stage 3 walkthrough](stage-3.md).
