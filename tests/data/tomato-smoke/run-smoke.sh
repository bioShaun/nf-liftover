#!/usr/bin/env bash
# Run the tomato 10kb smoke fixture.
#
# Usage:
#   bash tests/data/tomato-smoke/run-smoke.sh        # ID mode (default)
#   bash tests/data/tomato-smoke/run-smoke.sh id     # ID mode, explicit
#   bash tests/data/tomato-smoke/run-smoke.sh vcf    # VCF liftover mode
set -euo pipefail

MAMBA_BIN="${MAMBA_BIN:-/project/software/miniforge3/bin/mamba}"
NEXTFLOW_ENV="${NEXTFLOW_ENV:-nextflow26}"
NEXTFLOW_BIN="${NEXTFLOW_BIN:-nextflow}"

MODE="${1:-id}"
case "${MODE}" in
  id|vcf) ;;
  *) echo "Usage: $0 [id|vcf]" >&2; exit 2 ;;
esac

common_args=(
  run .
  --ref_fa tests/data/tomato-smoke/SL4.0ch01.10kb.fa
  --query_fa tests/data/tomato-smoke/LA2093.chr01.10kb.fa
  --mapping tests/data/tomato-smoke/chrom_pairs.tsv
  --outdir "results/tomato-smoke-${MODE}"
  -profile standard
  -resume
)

if [ "${MODE}" = "vcf" ]; then
  run_args=(
    "${common_args[@]}"
    --vcf tests/data/tomato-smoke/input.vcf
  )
else
  run_args=(
    "${common_args[@]}"
    --id tests/data/tomato-smoke/tomato-smoke.id
  )
fi

if command -v "${NEXTFLOW_BIN}" >/dev/null 2>&1; then
  exec "${NEXTFLOW_BIN}" "${run_args[@]}"
fi

exec "${MAMBA_BIN}" run -n "${NEXTFLOW_ENV}" nextflow "${run_args[@]}"
