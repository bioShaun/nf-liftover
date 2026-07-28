#!/usr/bin/env bash
# Run the tomato 10kb smoke fixture.
#
# Usage:
#   bash tests/data/tomato-smoke/run-smoke.sh        # ID mode (default)
#   bash tests/data/tomato-smoke/run-smoke.sh id     # ID mode, explicit
#   bash tests/data/tomato-smoke/run-smoke.sh vcf    # VCF liftover mode
#   bash tests/data/tomato-smoke/run-smoke.sh chain  # ID mode reusing existing chain
set -euo pipefail

MAMBA_BIN="${MAMBA_BIN:-/project/software/miniforge3/bin/mamba}"
NEXTFLOW_ENV="${NEXTFLOW_ENV:-nextflow26}"
NEXTFLOW_BIN="${NEXTFLOW_BIN:-nextflow}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "${ROOT}"

# Prefer committed example config when local nextflow.config is absent.
CONFIG_ARGS=()
if [ ! -f nextflow.config ] && [ -f nextflow.config.example ]; then
  CONFIG_ARGS=(-c nextflow.config.example)
fi

MODE="${1:-id}"
case "${MODE}" in
  id|vcf|chain) ;;
  *) echo "Usage: $0 [id|vcf|chain]" >&2; exit 2 ;;
esac

run_nf() {
  if command -v "${NEXTFLOW_BIN}" >/dev/null 2>&1; then
    "${NEXTFLOW_BIN}" "${CONFIG_ARGS[@]}" "$@"
  else
    "${MAMBA_BIN}" run -n "${NEXTFLOW_ENV}" nextflow "${CONFIG_ARGS[@]}" "$@"
  fi
}

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
elif [ "${MODE}" = "chain" ]; then
  CHAIN_SRC="${CHAIN_SRC:-results/tomato-smoke-id/chain/all.chain}"
  CHAIN_META_SRC="${CHAIN_META_SRC:-results/tomato-smoke-id/chain/chain_meta.yml}"
  if [ ! -f "${CHAIN_SRC}" ]; then
    echo "Missing chain at ${CHAIN_SRC}; run: $0 id" >&2
    exit 1
  fi
  run_args=(
    "${common_args[@]}"
    --id tests/data/tomato-smoke/tomato-smoke.id
    --chain "${CHAIN_SRC}"
  )
  if [ -f "${CHAIN_META_SRC}" ]; then
    run_args+=(--chain_meta "${CHAIN_META_SRC}")
  fi
else
  run_args=(
    "${common_args[@]}"
    --id tests/data/tomato-smoke/tomato-smoke.id
  )
fi

run_nf "${run_args[@]}"

OUT="results/tomato-smoke-${MODE}"
test -f "${OUT}/software_versions.yml"
test -f "${OUT}/run_meta.yml"
test -f "${OUT}/chain/chain_meta.yml"
test -f "${OUT}/chain/all.chain"
test -f "${OUT}/pipeline_info/execution_report.html"
! grep -qx 'END' "${OUT}/software_versions.yml"

# Strict path checks when shared env is used
if grep -E '_path:' "${OUT}/software_versions.yml" >/dev/null 2>&1; then
  bad=0
  while read -r key path; do
    case "${path}" in
      */envs/nf-liftover-tools/bin/*|*/conda/*|/envs/*/bin/*)
        echo "OK ${key}=${path}"
        ;;
      *)
        echo "BAD ${key}=${path}" >&2
        bad=1
        ;;
    esac
  done < <(grep -E '_path:' "${OUT}/software_versions.yml" | awk '{print $1, $2}' | sed 's/://')
  if [ "${bad}" -ne 0 ]; then
    exit 1
  fi
fi

if [ "${MODE}" = "chain" ] || [ "${MODE}" = "vcf" ]; then
  TRACE="${OUT}/pipeline_info/execution_trace.txt"
  if [ -f "${TRACE}" ]; then
    if grep -E 'ALIGN_WHOLE_CHROMOSOME|ALIGN_SPLIT_WINDOW|PAF_TO_CHAIN|COMBINE_ALL_PAFS' "${TRACE}"; then
      echo "alignment processes unexpectedly present in ${TRACE}" >&2
      exit 1
    fi
  fi
fi

echo "Smoke ${MODE} finished: ${OUT}"
