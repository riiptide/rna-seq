# Stage 2: executed notebook and cluster review

Completed November 4, 2026. `notebooks/01_explore.ipynb` ran from a fresh kernel
with all **14 code cells executed and zero errors**. The final run took 30.91
seconds with inputs/model cached (setup and download time excluded). A second
clean run reproduced the first run's Leiden and CellTypist labels exactly.

The notebook saved and reopened `results/annotated.h5ad`: **1,107 cells × 17,581
genes**, **13 Leiden clusters**, and **14 CellTypist label categories**, including
45 cells labeled Heterogeneous. The main output uses model majority-voted labels;
the working interpretations below are a separate human review, not silent overrides.

## How the analysis was chosen

- The barcode-rank curve drops steeply around 1,200 candidate barcodes. QC uses
  200–5,000 genes, 1,000–30,000 UMIs, and mitochondrial UMI fraction ≤20%.
  Candidate cells have a median mitochondrial fraction around 10%; a generic
  5% cutoff would exclude the main population. See the notebook's sensitivity table.
- The exact local Ensembl 108 mapping transforms 39,546 gene IDs into 39,509
  unique features. Thirty-six shared-symbol groups are summed before normalization.
  Missing symbols retain their IDs, and the full mapping table is exported.
- The main analysis retains genes detected in ≥3 cells, normalizes to 10,000 UMIs
  per cell, and applies log1p. Raw counts remain in `layers['counts']`.
- 1,998 nonmitochondrial HVGs feed scaled PCA. The strongest elbow is near PCs
  5–6; 20 PCs are retained conservatively to include subtler immune states. Those
  PCs explain 24.84% of scaled-HVG variance. The graph uses 15 neighbors and Leiden
  resolution 0.6, with seed 42. Cluster separation is not proof of distinct types.
- Immune_All_Low v2 receives a separate full-feature log1p(CP10k) matrix, not the
  scaled/HVG matrix. 5,842 of the model's 6,639 features overlap (88.0%). Model URL,
  observed SHA-256, parameters, and software versions are recorded.
- Majority voting uses graph overclusters at resolution 5.0 and minimum support
  0.5. It changes 155 labels, including 45 heterogeneous assignments. Both original
  predictions and voted labels are retained for comparison.

## What the clusters appear to be

These are marker-supported interpretations of this sample, not independently
validated identities. Fine subtypes and mixed or rare groups require caution.

| Cluster | Cells | Working interpretation | Evidence / caution |
|---|---:|---|---|
| 0 | 115 | Classical CD14+ monocytes, S100A8/A9-high | CD14, LYZ, S100A8, S100A9; 100% classical-monocyte model labels |
| 1 | 126 | Naive B cells | MS4A1, CD79A, TCL1A, IGHD/IGHM |
| 2 | 76 | Cytotoxic T cells; mixed effector/memory-like states | CD3E, CD8A/B, CCL5, NKG7; only 54% share the leading fine label |
| 3 | 174 | IL7R+ helper T cells; mixed naive/memory/effector-like states | CD3E, IL7R, TRAC; 63% leading fine label and 31 heterogeneous cells |
| 4 | 57 | Memory B cells | MS4A1, CD79A, TNFRSF13B, IGHG genes |
| 5 | 208 | Classical CD14+ monocytes, CPVL/LGALS2-high | CD14, LYZ, CPVL; 204 classical monocytes and 4 DC2 labels |
| 6 | 38 | Naive/central-memory-like CD8 T cells | CD3E, CD8A/B, CCR7, TCF7, LEF1 |
| 7 | 68 | MAIT-like T cells | CD3E, IL7R, NKG7, SLC4A10, KLRB1; all 68 raw predictions were MAIT |
| 8 | 57 | NK-like cytotoxic cells; review T-cell-associated RNA | GNLY, NKG7, KLRD1, FCER1G, TYROBP; CD3E/TRDC signal merits follow-up |
| 9 | 32 | Nonclassical FCGR3A+ monocytes | FCGR3A, CST3, LYZ, CDKN1C; weaker CD14 |
| 10 | 8 | Plasmacytoid dendritic cells (pDC), small group | 8 cells; LILRA4, CLEC4C, TCF4, GZMB; all raw predictions pDC |
| 11 | 130 | Naive/central-memory-like helper T cells | CD3E, IL7R, CCR7, TCF7, LEF1 |
| 12 | 18 | Conventional dendritic cells, DC2-like; subtype mixed | 18 cells; FCER1A, CST3, CD1C; raw labels include 11 DC2 and 5 DC1 |


Clusters 0 and 5 share the classical-monocyte lineage. Clusters 2 and 3 have
mixed T-cell model labels, so broad interpretations are more defensible than a
single confident subtype. Cluster 8 has strong NK-associated genes and little
CD3D/CD8B, but its CD3E/TRDC signal warrants follow-up. TRDC alone cannot establish
a gamma-delta T-cell identity. Cluster 12's majority vote assigns DC2 throughout,
while raw predictions include DC1 and other dendritic labels. The eight-cell pDC
cluster is marker-supported but too small to claim robust population estimates.

The requested twelve PBMC markers were all present. Additional markers in the
notebook help assess the MAIT-like, pDC, naive-like T-cell, and NK-like groups.
Reference context: [Azimuth human PBMC markers](https://azimuth.hubmapconsortium.org/references/human_pbmc/)
and [Human Protein Atlas TRDC](https://www.proteinatlas.org/ENSG00000211829-TRDC/single%2Bcell).

## Artifacts

- Notebook: `notebooks/01_explore.ipynb`, with executed figures and tables.
- Final AnnData: `results/annotated.h5ad`.
- QC decisions and mapping audit: `results/stage2/barcode_qc.csv.gz` and
  `results/stage2/gene_id_mapping.tsv`.
- Ranked genes: `results/stage2/marker_genes_all.csv` and `marker_genes_top.csv`.
- Cluster evidence: `results/stage2/cluster_summary.csv`, before/after voting
  crosstabs, and marker expression/fraction tables.
- Figures: `results/stage2/figures/`.
- [Execution record](stage-2-execution.json), [run provenance](stage-2-run.json),
  and [review table](stage-2-cluster-review.csv).
- [Environment and launch instructions](stage-2.md).

The exported AnnData was checked for valid sparse counts, matching expression
values and identifiers after reopening, finite UMAP coordinates, labels for every
cell, a neighbor graph, and marker rankings. Original Stage 1 inputs are unchanged.
Two informational dependency warnings remain: scaling densifies only the small
HVG copy, and CellTypist accesses a deprecated Scanpy version property. Neither
interrupted execution. Plot layouts were inspected.

No formal empty-droplet test, dedicated doublet detection, or ambient-RNA
correction was performed. QC ceilings do not establish a singlet population.
One donor and cluster-derived marker tests do not support population-level
condition comparisons. Model scores are not calibrated probabilities; the score
panel is the original cell-level best score, not confidence in a changed voted label.

Stage 2 stops here for review. Nextflow modules and Docker are Stage 3.
