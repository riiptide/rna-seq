# Stage 4: explore the saved analysis

The dashboard turns the Stage 3 output into an interactive view. It reads saved
results without rerunning clustering, annotation, or filtering.

## Setup

Use Python 3.11 and run these commands from the repository root. A separate
virtual environment keeps dashboard dependencies independent of the analysis
environments and Docker images.

```bash
python3.11 -m venv work/envs/dashboard
work/envs/dashboard/bin/python -m pip install -r dashboard/requirements.lock.txt
work/envs/dashboard/bin/streamlit run dashboard/app.py
```

Open <http://localhost:8501>. Stop the server with Ctrl+C. The default configuration
binds to localhost; the phone layout can be previewed using a narrow browser window.

The default inputs are `results/analysis/annotated.h5ad` and
`results/analysis/markers.csv`, produced by the full pipeline. To explore another
run, supply both matching files after Streamlit's `--` separator:

```bash
work/envs/dashboard/bin/streamlit run dashboard/app.py -- \
  --data results/test-pipeline/analysis/annotated.h5ad \
  --markers results/test-pipeline/analysis/markers.csv
```

If the matrix is missing, the app explains how to generate it. A missing marker
table produces a warning in that tab while the other views remain available.
Read-depth subsampling makes test-profile labels unsuitable for biological validation.

## What each piece does

- **Streamlit** builds the controls, tabs, summary cards, and table from Python.
  Its caching avoids rereading the matrix when a control changes; file size and
  modification time invalidate the cache when a pipeline output is replaced.
- **Plotly** draws interactive UMAPs and histograms with hover, pan, zoom, reset,
  and PNG export. Stable category colors make repeated comparisons easier.
- **AnnData** reads the `.h5ad`, which holds cell metadata, expression, the UMAP
  coordinates, and the analysis parameters together.
- **pandas, NumPy, and SciPy** select metadata and expression values. Only the
  selected gene's sparse expression vector is converted to a dense array.
- **CSS** controls presentation. Below 720 pixels, cards, plots, and legends stack
  vertically; the marker table can scroll horizontally within its panel.

`dashboard/app.py` connects the controls to data and charts. `dashboard/data.py`
validates and extracts saved values, and `dashboard/charts.py` constructs the
figures. `dashboard/style.css` and `.streamlit/config.toml` define presentation.
`requirements.txt` pins direct dependencies; `requirements.lock.txt` also pins
transitive dependencies from the tested Python 3.11 environment on macOS arm64.
Other platforms have not been verified.

## Explore the sample

1. Start with **Cell type** to see CellTypist Immune_All_Low majority-voted labels.
   The full result contains **1,107 cells**, **13 Leiden clusters**, a median of
   **2,066 genes per cell**, and **13 predicted cell types**. Another **45 cells**
   are labeled Heterogeneous and remain included in the cell total.
2. Switch to **Leiden cluster** to compare expression-based groups with predicted
   identities. A cluster and a cell type need not have a one-to-one relationship.
3. Choose **Gene expression** and type a gene symbol. Try `MS4A1` for B cells and
   `CD3E` for T cells. Colors use `.raw.X`, which this pipeline stores after
   normalization to 10,000 UMIs per cell and `log1p`; this is not raw counts.
   The plot draws cells with higher expression last so they remain visible.
   MS4A1 is detected in 200 cells (18.1%); CD3E in 514 cells (46.4%).
4. Open **Quality control**. Histograms show retained cells, not the rejected
   population. Dashed lines come from the saved parameters: 200–5,000 genes,
   1,000–30,000 UMIs, and at most 20% mitochondrial UMIs. Nine retained cells have
   a Scrublet doublet flag; Stage 3 flags doublets but does not remove them by default.
   Expand the settings to inspect the run's recorded choices.
5. Open **Marker genes**, select a cluster, and compare its leading predicted
   label with enriched genes. Positive Wilcoxon scores are sorted descending;
   adjusted p-values come from Scanpy's Benjamini–Hochberg correction. The row
   limit controls the displayed table. **Download cluster markers** exports all
   positive-score rows for the selected cluster, including the original CSV fields.

The dashboard reports the pipeline's annotations; it does not independently
validate cell identities. Gene names without a reference symbol retain their
Ensembl identifiers. It expects the Stage 3 schema, including log-normalized
expression in `.raw`, and a marker CSV from the same run.

## Stage 4 verification

Verified against the full Stage 3 output on November 23, 2026:

- Python imports/compilation and `pip check` passed in the isolated environment.
- Streamlit AppTest loaded the real matrix and marker table without exceptions;
  the four metrics matched the saved data. Gene expression, gene selection,
  Leiden coloring, and cluster selection were exercised. Missing matrix and marker
  files produced the expected messages without uncaught exceptions.
- Source checks confirmed gene detection counts, histogram cell totals, and
  matching cluster identifiers between the matrix and marker table.
- Browser checks confirmed cell-type and Leiden coloring, a CD3E gene search,
  recorded QC lines, and the cluster-marker table.
- Layouts were inspected at 1360 × 1000 and 390 × 844. Phone-width cards and QC
  plots stack, controls wrap, and the page has no horizontal overflow. This is
  browser viewport verification, not testing on a physical phone.

These are Stage 4 smoke checks. The planned pytest suite and GitHub Actions
remain Stage 5.
