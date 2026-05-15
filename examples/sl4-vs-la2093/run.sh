#!/usr/bin/env bash
set -euo pipefail

NEXTFLOW_BIN="${NEXTFLOW_BIN:-nextflow}"
if ! command -v "${NEXTFLOW_BIN}" >/dev/null 2>&1; then
  NEXTFLOW_BIN="/project/software/miniforge3/envs/nextflow/bin/nextflow"
fi

"${NEXTFLOW_BIN}" run /public/scripts/nf-liftover \
  --id        /public/data/genomes/solanum_lycopersicum_LA2093/TCZZSL20K.id \
  --ref_fa    /public/data/genomes/solanum_lycopersicum_LA2093/S_lycopersicum_chromosomes.4.00.fa \
  --query_fa  /public/data/genomes/solanum_lycopersicum_LA2093/rename.genome.fa \
  --mapping   /public/scripts/nf-liftover/examples/sl4-vs-la2093/chrom_pairs.tsv \
  --outdir    results/TCZZSL20K-LA2093 \
  -profile    standard \
  -resume
