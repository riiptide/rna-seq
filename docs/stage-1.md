# Stage 1: from reads to counts

Run these commands from the repository root, in **Bash** with Conda initialized.
This walkthrough is specific to the verified Apple Silicon Mac setup (24 GB RAM).
The November 4 setup reused checksum-verified downloads from an earlier local
Stage 1 project and recreated environments from the locks below. The download
commands remain the instructions for reproducing the project on a new machine.
It invokes FastQC and kb directly; Nextflow only reads configuration at this stage.

## 1. Understand the inputs

Each FASTQ record has four lines: read identifier, sequence, separator, and base
quality scores. In 10x 3′ v3, R1 contains a 16-base cell barcode and a 12-base UMI
(unique molecular identifier). R2 contains the transcript sequence. I1 holds a
sample index; these public reads are already demultiplexed, so counting uses R1/R2.

A barcode identifies a droplet. A UMI helps distinguish original RNA molecules
from PCR copies. A barcode with reads is not automatically a real cell. Both
sequencing lanes belong to the same sample; count them together to deduplicate
molecules across lanes.

## 2. Create the isolated environments

Miniforge/Mamba installs software into an environment without changing your base
Python installation. This setup deliberately pins compatible versions, including
Python 3.11.11, kb-python 0.30.2, pandas 2.2.3, AnnData 0.11.4, FastQC 0.12.1,
seqtk 1.4, Java 17.0.18, and bustools 0.45.1.

```bash
export MAMBA_ROOT_PREFIX="$PWD/work/mamba"
export CONDA_PKGS_DIRS="$PWD/work/cache/conda"
export PIP_CACHE_DIR="$PWD/work/cache/pip"
mamba env create -p "$PWD/work/envs/stage1" -f environment.yml
mamba env create -p "$PWD/work/envs/native-count" -f environment-native.yml
conda activate "$PWD/work/envs/stage1"
python -m pip check
```

For an exact repeat on **osx-arm64**, use these commands instead of the two
`mamba env create` commands above:

```bash
mamba create -y -p "$PWD/work/envs/stage1" --file docs/conda-osx-arm64.lock.txt
work/envs/stage1/bin/python -m pip install -r docs/requirements-stage1.lock.txt
mamba create -y -p "$PWD/work/envs/native-count" --file docs/conda-native-osx-arm64.lock.txt
```

Do not create a second environment over an existing one. Activate the supplied
local environment to inspect this run; use a fresh clone for a full reproduction.
The lock files record exact packages; they do not promise cross-platform binaries.
Linux containers will be introduced and verified in Stage 3.

### The Mac compatibility detail

The working kallisto executable is the **standard-k Apple Silicon binary bundled
with kb-python**, version 0.51.1. Its HDF5 dependency is supplied by the second
local environment (HDF5 1.14.3). The count commands below explicitly supply its
library search path. They also explicitly select bustools 0.45.1.

The tested Bioconda kallisto 0.51.1 and 0.52.0 builds aborted when reading this
prebuilt index. The current Bioconda recipe enables `MAX_KMER_SIZE=64`; build
options can matter as well as version numbers. The bundled standard-k binary
was verified with `kallisto inspect` and actual reads. We retain the downloaded
prebuilt reference, rather than rebuilding it silently.

Pandas 3.1.0 also caused AnnData's nullable-string export guard to reject the
barcode metadata. Pandas 2.2.3 is pinned to keep the verified export behavior.
These were environment compatibility issues, not evidence of bad biological data.

## 3. Load configuration

Nextflow 26.10.5 is used here as a configuration reader, without running processes.
Java is supplied by the activated environment. The launcher is already installed
on the development Mac. For another machine, follow the official Nextflow
installation guide linked below and select the same version.

```bash
export NXF_VER=26.10.5
export NXF_HOME="$PWD/work/nextflow"
export MPLCONFIGDIR="$PWD/work/cache/matplotlib"
export NUMBA_CACHE_DIR="$PWD/work/cache/numba"
source docs/load-params.sh
```

`load-params.sh` asks Nextflow to export `nextflow.config` as JSON, then uses
Python's `shlex.quote` to turn scalar values into shell variables safely. For
example, `params.seed` becomes `$SEED`. Paths, chemistry, threads, and sampling
settings are defined once. Rerun this command after changing configuration.

## 4. Download, verify, and extract

```bash
mkdir -p work/downloads work/logs "$RAW_DIR" "$TEST_DIR" "$REFERENCE_DIR"
curl --fail --location --retry "$DOWNLOAD_RETRIES" --continue-at - \
  --output "$RAW_ARCHIVE" "$RAW_URL"
curl --fail --location --retry "$DOWNLOAD_RETRIES" --continue-at - \
  --output "$REFERENCE_ARCHIVE" "$REFERENCE_URL"
shasum -a 256 -c docs/downloads.sha256
tar -tf "$RAW_ARCHIVE"
tar -tf "$REFERENCE_ARCHIVE"
tar -xf "$RAW_ARCHIVE" -C "$RAW_DIR"
tar -xf "$REFERENCE_ARCHIVE" -C "$REFERENCE_DIR"
shasum -a 256 -c docs/inputs.sha256
```

`curl --fail` detects HTTP errors, `--location` follows redirects, `--retry` handles
transient failures, and `--continue-at -` resumes partial downloads. The FASTQ
archive is 5,549,312,000 bytes; the reference archive is 145,139,492 bytes.

`shasum` checks that bytes match this recorded run. These are locally recorded
SHA-256 values, not independent checksums published by 10x. `tar -tf` lists an
archive before extraction. These verified archives contain the relative paths
shown in `samplesheet.csv` and the reference directory.

The reference is human GRCh38 / Ensembl 108, standard transcriptome index,
Pachter lab release v1. `index.idx` supports fast transcript compatibility lookup;
`t2g.txt` maps transcript identifiers to gene identifiers and symbols. They must
come from the same reference release. The standard reference includes mature
transcripts from the biotypes documented by the reference publisher.

## 5. Make and validate the test reads

```bash
set -euo pipefail
while IFS=, read -r sample lane r1 r2; do
  seqtk sample -s"$SEED" "$r1" "$TEST_PAIRS_PER_LANE" |
    gzip -n > "$TEST_DIR/${sample}_${lane}_R1.fastq.gz"
  seqtk sample -s"$SEED" "$r2" "$TEST_PAIRS_PER_LANE" |
    gzip -n > "$TEST_DIR/${sample}_${lane}_R2.fastq.gz"
  python bin/validate_fastq_pairs.py \
    "$TEST_DIR/${sample}_${lane}_R1.fastq.gz" \
    "$TEST_DIR/${sample}_${lane}_R2.fastq.gz" \
    --expected-pairs "$TEST_PAIRS_PER_LANE"
done < <(tail -n +2 "$SAMPLESHEET")
shasum -a 256 -c docs/test-inputs.sha256
```

The seed is 42 and the initial fixture has 500,000 pairs per lane (1,000,000 total).
The two input mates must already have the same read order; the same seed makes
seqtk choose matching positions. The validator checks every sampled record's
structure and identifier, and the exact pair count. `gzip -n` omits the original
filename and timestamp from its header. Exact compressed bytes can still depend
on compressor version; the recorded checksums describe this Mac's output.

`set -euo pipefail` makes failures visible, including failures upstream in a pipe.
The `while` loop reads the simple, unquoted CSV used by this repository. Its
columns are sample, lane, R1 path, R2 path. A general quoted-CSV parser belongs in
the future workflow; don't insert commas into these fields.

Subsampling reduces depth, not the number of cells. This fixture tests execution;
its sparse expression profiles must not be used to claim annotation accuracy.
Whether it supports the full downstream test in under five minutes will be
measured in Stage 3. Fixture preparation scans the full FASTQs and is excluded
from that timed execution.

## 6. Assemble the two lane pairs

```bash
FULL_READS=()
TEST_READS=()
while IFS=, read -r sample lane r1 r2; do
  FULL_READS+=("$r1" "$r2")
  TEST_READS+=("$TEST_DIR/${sample}_${lane}_R1.fastq.gz"
               "$TEST_DIR/${sample}_${lane}_R2.fastq.gz")
done < <(tail -n +2 "$SAMPLESHEET")
```

Bash arrays keep each path as a separate argument. Ordering is `L001_R1 L001_R2
L002_R1 L002_R2`; each R2 immediately follows its R1. The samplesheet currently
contains one biological sample. Future Nextflow modules will group lanes by sample.

## 7. Run FastQC manually

```bash
mkdir -p "$FULL_FASTQC" "$TEST_FASTQC"
/usr/bin/time -l fastqc --threads "$FASTQC_THREADS" \
  --memory "$FASTQC_MEMORY_MB" --outdir "$TEST_FASTQC" "${TEST_READS[@]}" \
  > work/logs/fastqc-test.log 2> work/logs/fastqc-test.time
/usr/bin/time -l fastqc --threads "$FASTQC_THREADS" \
  --memory "$FASTQC_MEMORY_MB" --outdir "$FULL_FASTQC" "${FULL_READS[@]}" \
  > work/logs/fastqc-full.log 2> work/logs/fastqc-full.time
```

FastQC reports base qualities, sequence composition, adapter content, duplication,
and other diagnostics. It does **not** trim or alter the reads. It writes an HTML
report for inspection and a ZIP containing machine-readable results. The configured
512 MB allocation is per worker; four input files may be processed together.

Open the HTML files under `results/full/fastqc/` in your browser. The compact
[status table](fastqc-summary.tsv) includes all eight reports from this run.
Base-quality and adapter checks pass. Composition warnings occur on all files;
full-depth files also flag sequence duplication. Single-cell barcode structure,
RNA abundance, and PCR duplication affect these diagnostics. These flags alone
are not a reason to trim R1 or remove duplicate reads before UMI-aware counting.

## 8. Run kb count manually

Run the test count first, validate it, then run the full count:

```bash
/usr/bin/time -l env DYLD_FALLBACK_LIBRARY_PATH="$PWD/$COUNT_LIBRARY_PATH" \
  kb count --kallisto "$KALLISTO" --bustools "$BUSTOOLS" \
  -i "$INDEX" -g "$T2G" -x "$CHEMISTRY" --workflow "$WORKFLOW" \
  -t "$THREADS" -m "$SORT_MEMORY" --h5ad \
  -o "$TEST_COUNTS" "${TEST_READS[@]}" \
  > work/logs/kb-test.log 2> work/logs/kb-test.time
python bin/inspect_counts.py "$TEST_COUNTS/counts_unfiltered/adata.h5ad"

/usr/bin/time -l env DYLD_FALLBACK_LIBRARY_PATH="$PWD/$COUNT_LIBRARY_PATH" \
  kb count --kallisto "$KALLISTO" --bustools "$BUSTOOLS" \
  -i "$INDEX" -g "$T2G" -x "$CHEMISTRY" --workflow "$WORKFLOW" \
  -t "$THREADS" -m "$SORT_MEMORY" --h5ad \
  -o "$FULL_COUNTS" "${FULL_READS[@]}" \
  > work/logs/kb-full.log 2> work/logs/kb-full.time
python bin/inspect_counts.py "$FULL_COUNTS/counts_unfiltered/adata.h5ad"
```

The `env` assignment after `/usr/bin/time` makes the Mac HDF5 library path available
to the counting process. It applies to this command, not the system environment.
`/usr/bin/time -l` records elapsed time and process resource statistics on macOS;
Linux uses a different `time` interface. This timing is a baseline, not the optional
benchmark add-on. The recorded maximum RSS is not aggregate simultaneous memory
across every process in a future workflow.

| Argument | Meaning |
|---|---|
| `--kallisto`, `--bustools` | Select the exact binaries that were verified |
| `-i`, `-g` | Matching transcriptome index and transcript-to-gene mapping |
| `-x 10xv3` | Decode the v3 barcode/UMI/read structure |
| `--workflow standard` | Quantify mature transcript signal |
| `-t` | Number of compute threads |
| `-m 4G` | BUS sorting-memory allocation, not a total-job memory limit |
| `--h5ad` | Export the sparse raw count matrix as an AnnData file |
| `-o` | Output directory |

Kallisto pseudoaligns transcript reads using compatible sequence k-mers. Its BUS
output associates **B**arcodes, **U**MIs, and transcript compatibility **S**ets.
Bustools sorts records, corrects barcodes against the packaged 10x v3 on-list,
groups compatible transcripts into genes, and collapses UMI duplicates.

The on-list describes valid assay barcodes; it is not a list of confirmed cells.
We omit `--filter` to preserve the unfiltered matrix for examining cell calling
in Stage 2. Gene identifiers are retained; `t2g.txt` supplies gene symbols for
mitochondrial QC, marker interpretation, and later CellTypist annotation.
The count log warns about 13,914 genes with blank or duplicated symbols in the
reference. Their gene IDs remain intact. Stage 2 will explicitly handle these
symbol cases before annotation; we do not silently remove these genes here.

**Validate the output, not just the exit status.** In this version, kb can log an
exception yet return zero. A partially written `.h5ad` can also exist after an
export failure. `inspect_counts.py` reopens it and checks a nonempty sparse matrix,
unique identifiers, finite nonnegative integer-valued counts, and nonzero entries.
Stage 3 must preserve these checks around the counting process.

For a fresh timed rerun, move the previous count directory aside first. kb can
reuse existing intermediate files; measuring that would not be a clean count run.

## 9. Inspect and review

Expected outputs inside each count directory:

- `run_info.json`: reads processed and pseudoalignment statistics.
- `inspect.json`: BUS-level barcode/read diagnostics before correction.
- `counts_unfiltered/cells_x_genes.mtx`: sparse raw UMI counts.
- Matching `.barcodes.txt` and `.genes.txt` files: matrix row/column identifiers.
- `counts_unfiltered/adata.h5ad`: the same raw counts packaged as AnnData.

See [the run record](stage-1-results.md) for measured results and
[tool versions](tool-versions.json) for the final environment. Full reads,
references, environments, and generated matrices are ignored by Git. Small
provenance records and configuration are committed.

Before Stage 2, review these questions:

1. Why does trimming R1 risk destroying the barcode/UMI layout?
2. Why must lanes from one sample be counted together?
3. Why does an unfiltered matrix have far more rows than the vendor's 1,222 cells?
4. Why can shallow subsampling change the cells that pass QC?

Stop here. Normalization, cell QC, clustering, UMAP, marker testing, and annotation
belong to the next stage. The `.h5ad` here is an intermediate raw-count artifact,
not the final annotated analysis output.

## Sources

- [10x dataset and CC BY 4.0 attribution](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0)
- [Reference release and construction](https://github.com/pachterlab/kallisto-transcriptome-indices)
- [kb-python manual](https://kallisto.readthedocs.io/en/latest/manuals/kb.html)
- [seqtk paired subsampling](https://github.com/lh3/seqtk)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Nextflow installation](https://docs.seqera.io/nextflow/install)
- [Bioconda kallisto build recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/kallisto/build.sh)
