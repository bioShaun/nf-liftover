#!/usr/bin/env bash
set -euo pipefail

NEXTFLOW_BIN="${NEXTFLOW_BIN:-nextflow}"
if ! command -v "${NEXTFLOW_BIN}" >/dev/null 2>&1; then
  NEXTFLOW_BIN="/project/software/miniforge3/envs/nextflow26/bin/nextflow"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 注意（2026-08-05 审计，见 docs/drift-audit-2026-08-05.md）：
# 下方 --id / --ref_fa / --query_fa 引用的目录
#   /public/data/genomes/solanum_lycopersicum_LA2093/
# 在当前主机上不存在（复数和单数 genome 路径均已核实缺失）。
# 运行前请：
#   1) 自备 --id / --ref_fa / --query_fa 指向实际可用的文件；或
#   2) 恢复 LA2093 数据集到上述路径。
# --mapping 指向的 chrom_pairs.tsv 在部署副本中存在，可保留。
# ─────────────────────────────────────────────────────────────────────────────

"${NEXTFLOW_BIN}" run /public/scripts/nf-liftover \
  --id        /public/data/genomes/solanum_lycopersicum_LA2093/TCZZSL20K.id \
  --ref_fa    /public/data/genomes/solanum_lycopersicum_LA2093/S_lycopersicum_chromosomes.4.00.fa \
  --query_fa  /public/data/genomes/solanum_lycopersicum_LA2093/rename.genome.fa \
  --mapping   /public/scripts/nf-liftover/examples/sl4-vs-la2093/chrom_pairs.tsv \
  --outdir    results/TCZZSL20K-LA2093 \
  -profile    standard \
  -resume
