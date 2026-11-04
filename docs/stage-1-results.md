# Stage 1 review: reads to counts

Completed November 4, 2026 (America/New_York), in the Desktop `rna-seq` repo,
on an Apple Silicon Mac with 24 GiB RAM (verified using `sysctl`).

## What was done

- Initialized Git and connected `origin` to `https://github.com/riiptide/rna-seq.git`.
- Reused an earlier local Stage 1 implementation (source commit `6a0ef72`) and
  input downloads; verified archive and extracted-input SHA-256 checksums.
- Recreated both isolated environments from their package locks. Fixed a
  nonportable build-directory URL in the pip lock (`packaging==26.3`).
- Verified pinned package versions and kallisto/bustools executable checksums.
- Regenerated the test reads with seqtk and checked every sampled mate pair.
- Ran FastQC and kb count directly from Bash into fresh output directories.
- Reopened both AnnData matrices and checked the data, not just command exit codes.

The [walkthrough](stage-1.md) explains the commands. The [provenance record](provenance.json)
distinguishes reused downloads from this fresh run. Checksums are locally recorded;
they are not independent checksums published by the data provider.

## Measured results

| Measurement | Test subset | Full sample |
|---|---:|---:|
| Input read pairs | 1,000,000 | 66,601,887 |
| Pseudoaligned read pairs | 561,276 | 37,366,584 |
| Pseudoalignment (%) | 56.13 | 56.10 |
| FastQC elapsed (seconds) | 6.62 | 174.24 |
| kb count elapsed (seconds) | 41.83 | 396.77 |
| Count maximum RSS (GiB) | 4.02 | 4.03 |
| Unfiltered barcodes | 23,175 | 237,012 |
| Reference genes | 39,546 | 39,546 |
| Genes with counts | 13,453 | 22,577 |
| Counted UMIs | 503,795 | 10,920,109 |
| Zero-count barcodes | 955 | 4,548 |

These rows are **unfiltered barcodes, not confirmed cells**. Empty droplets and
low-quality barcodes remain; selecting cells is Stage 2 work. Counts are sparse,
finite, nonnegative, integer-valued UMIs stored using a floating-point dtype.
Gene and barcode identifiers are unique. Zero-count rows will be removed during
cell QC, not silently treated as cells here.

## Files to inspect

- `results/full/kb/counts_unfiltered/adata.h5ad`: full raw counts for Stage 2.
- `results/test/kb/counts_unfiltered/adata.h5ad`: shallow test counts.
- `results/full/fastqc/*_fastqc.html`: full-depth read-quality reports.
- [FastQC module statuses](fastqc-summary.tsv).
- [Full count validation](counts-full.json) and [test count validation](counts-test.json).
- [Machine-readable timings and counts](run-metrics.json).

Generated matrices, reads, references, and environments are local and ignored by
Git. Small source files, provenance, locks, and measured run records are committed.
Commits remain local pending your review; the GitHub remote is configured.

## Interpretation and next-stage decisions

Base-quality and adapter checks passed in all eight reports. Sequence composition
was flagged in all files, and full-depth duplication checks failed for lane 1
and warned for lane 2. FastQC flags are diagnostics. R1 contains barcodes and UMIs, and RNA expression
and amplification produce nonuniform sequences. Do not trim the v3 R1 barcode
layout or discard duplicate reads based only on these flags. UMI-aware counting
handles molecule duplication. Inspect [all statuses](fastqc-summary.tsv) alongside
the plots rather than using an all-green report as the acceptance criterion.

The reference has 227,665 transcripts and 39,546 gene IDs. Some genes lack symbols
or share symbols: 13,841 gene IDs have blank symbols, and 73 IDs have duplicated
nonempty symbols. Stage 2 must preserve gene IDs and resolve the symbols needed
for CellTypist. This standard transcriptome reference excludes intronic signal
and some gene biotypes; mitochondrial percentages will refer to this reference.

The test fixture contains 500,000 read pairs per lane, sampled with seed 42.
It tests execution; its shallow expression profiles cannot establish biological
annotation accuracy. We will assess whether it retains sufficient signal for
the downstream test profile when implementing those steps.

## Timing limits

Test FastQC plus counting took **48.45 seconds** summed across the two commands.
This excludes downloads, environment setup, fixture preparation, Scanpy,
CellTypist, and Nextflow/container overhead. The full under-five-minute pipeline
target is not yet established. Full FastQC overlapped fixture preparation;
these are practical run measurements, not a controlled performance benchmark.
Maximum RSS comes from macOS `/usr/bin/time -l`, not aggregate memory across all
simultaneous workflow processes. The kb `-m 4G` option limits sorting allocation,
not total job memory.

## Verification and review

Eleven focused validator checks passed, covering correct reads/counts and rejecting
mate mismatches, malformed/empty reads, empty identifiers, fractional/negative/
nonfinite counts, all-zero matrices, and dense matrices. `pip check` found no
dependency conflicts. Matrix Market values and row/column identifiers match each
exported AnnData matrix exactly; output SHA-256 values are recorded. The formal
pytest suite and CI belong to Stage 5.

Before Stage 2, review:

1. What information does R1 carry, and why would trimming it be risky?
2. Why do we count both lanes together?
3. Why are there more barcode rows than the vendor's reported 1,222 cells?
4. Why does a tiny random read subset change which cells pass QC?

**Stop here for review.** The notebook, annotation, Nextflow modules, containers,
dashboard, and CI have not been implemented. The current `.h5ad` is an intermediate
raw-count artifact, not the final annotated output.
