#!/usr/bin/env bash
set -euo pipefail

MAMBA_BIN="${MAMBA_BIN:-/project/software/miniforge3/bin/mamba}"
NEXTFLOW_ENV="${NEXTFLOW_ENV:-nextflow26}"
NEXTFLOW_BIN="${NEXTFLOW_BIN:-nextflow}"

run_args=(
  run .
  --id tests/data/tomato-smoke/tomato-smoke.id
  --ref_fa tests/data/tomato-smoke/SL4.0ch01.10kb.fa
  --query_fa tests/data/tomato-smoke/LA2093.chr01.10kb.fa
  --mapping tests/data/tomato-smoke/chrom_pairs.tsv
  --outdir results/tomato-smoke
  -profile standard
  -resume
)

if command -v "${NEXTFLOW_BIN}" >/dev/null 2>&1; then
  exec "${NEXTFLOW_BIN}" "${run_args[@]}"
fi

exec "${MAMBA_BIN}" run -n "${NEXTFLOW_ENV}" nextflow "${run_args[@]}"
