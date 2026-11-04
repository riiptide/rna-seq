#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
docker build --platform linux/amd64 -t rna-seq-analysis:0.3.0 containers/analysis
docker build --platform linux/amd64 -t rna-seq-kb:0.30.2-k0.51.1-1 containers/kb
docker pull --platform linux/amd64 quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0@sha256:e194048df39c3145d9b4e0a14f4da20b59d59250465b6f2a9cb698445fd45900
docker pull --platform linux/amd64 quay.io/biocontainers/multiqc:1.27--pyhdfd78af_0@sha256:a024b23be16fa3f3919ef018741478a7a3bfaaed0222369dff59abaf28316149
