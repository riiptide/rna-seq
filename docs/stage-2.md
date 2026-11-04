# Stage 2: interactive analysis

Open `notebooks/01_explore.ipynb`. Every code cell has a preceding explanation.
This stage runs on the full Stage 1 counts; it does not implement Nextflow modules.

## Open the existing local environment

From the repository root:

```bash
export JUPYTER_PATH="$PWD/work/envs/stage2/share/jupyter"
export JUPYTER_RUNTIME_DIR="$PWD/work/jupyter/runtime"
export IPYTHONDIR="$PWD/work/ipython"
mkdir -p "$JUPYTER_RUNTIME_DIR" "$IPYTHONDIR"
work/envs/stage2/bin/jupyter notebook notebooks/01_explore.ipynb
```

Select **RNA-seq Stage 2 (Python 3.11)**. Use **Restart Kernel and Run All Cells**
to reproduce the analysis, or run cells individually while learning. The notebook
finds the repository root from its working directory. Output paths and analysis
parameters are read from `nextflow.config`; edit those there and rerun from the top.

## Recreate the environment on Apple Silicon

Use a new environment directory rather than overwriting an existing environment:

```bash
export MAMBA_ROOT_PREFIX="$PWD/work/mamba"
export CONDA_PKGS_DIRS="$PWD/work/cache/conda"
export PIP_CACHE_DIR="$PWD/work/cache/pip"
mamba create -y -p "$PWD/work/envs/stage2" --file docs/conda-osx-arm64.lock.txt
work/envs/stage2/bin/python -m pip install -r docs/requirements-stage2.lock.txt
work/envs/stage2/bin/python -m pip check
work/envs/stage2/bin/python -m ipykernel install --prefix "$PWD/work/envs/stage2" \
  --name rna-seq-stage2 --display-name 'RNA-seq Stage 2 (Python 3.11)'
```

The Conda lock records the tested osx-arm64 builds. `environment-stage2.yml` is the
human-readable environment specification; the full Python lock records the actual
resolved versions. Stage 1 retains its own environment. Nextflow 26.10.5 remains a
configuration reader here and uses the environment's Java runtime. Linux/Docker
portability will be verified in Stage 3.

Core versions: Python 3.11.11, Scanpy 1.11.5, AnnData 0.11.4, CellTypist 1.7.1,
igraph 0.11.9, leidenalg 0.10.2, ipykernel 6.29.5, notebook 7.3.3.

## Inputs and outputs

- Input: `results/full/kb/counts_unfiltered/adata.h5ad` and matching `t2g.txt`.
- Model: official Immune_All_Low v2; download URL and observed SHA-256 are pinned
  in `nextflow.config`. A missing model is downloaded; a mismatched model is rejected.
- Result: `results/annotated.h5ad`.
- Supporting tables and figures: `results/stage2/`.
- All-barcode QC decisions: `results/stage2/barcode_qc.csv.gz`.
- Marker tables: `results/stage2/marker_genes_all.csv` and `marker_genes_top.csv`.
- Per-cluster labels and evidence: `results/stage2/cluster_summary.csv` and
  `cluster_celltypist_crosstab.csv`.

`layers['counts']` contains raw UMIs; `.X` and `.raw.X` contain log-normalized
expression for analysis genes. The model receives a separate normalized copy with
all mapped features, including those removed by the analysis gene filter.
Unmapped symbols retain Ensembl IDs, duplicate symbols are summed before
normalization, and the mapping audit preserves original gene IDs.

## Interpretation

QC thresholds are exploratory choices informed by this sample's plots. A UMI or
gene ceiling does not establish that all doublets were removed. Mitochondrial
percentage is based on quantified UMI counts, not raw read counts. The standard
reference excludes introns and some gene biotypes. This is one donor: marker-test
p-values do not establish population-level differential expression.

CellTypist predicts fine labels and smooths them within small graph subclusters.
Its scores are not calibrated clinical probabilities. Inspect marker combinations,
label mixtures, and cluster sizes before accepting a label; changing parameters
can change cluster numbering. Review the executed notebook and `stage-2-results.md`
before starting Stage 3.
