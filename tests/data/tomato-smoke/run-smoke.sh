#!/usr/bin/env bash
set -euo pipefail

NEXTFLOW_BIN="${NEXTFLOW_BIN:-nextflow}"
if ! command -v "${NEXTFLOW_BIN}" >/dev/null 2>&1; then
  NEXTFLOW_BIN="/project/software/miniforge3/envs/nextflow/bin/nextflow"
fi

"${NEXTFLOW_BIN}" run . \
  --id tests/data/tomato-smoke/tomato-smoke.id \
  --ref_fa tests/data/tomato-smoke/SL4.0ch01.10kb.fa \
  --query_fa tests/data/tomato-smoke/LA2093.chr01.10kb.fa \
  --mapping tests/data/tomato-smoke/chrom_pairs.tsv \
  --outdir results/tomato-smoke \
  -profile standard \
  -resume
