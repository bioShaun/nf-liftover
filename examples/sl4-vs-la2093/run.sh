#!/usr/bin/env bash
set -euo pipefail

nextflow run /public/scripts/nf-liftover \
  --id        /public/data/genomes/solanum_lycopersicum_LA2093/TCZZSL20K.id \
  --ref_fa    /public/data/genomes/solanum_lycopersicum_LA2093/S_lycopersicum_chromosomes.4.00.fa \
  --query_fa  /public/data/genomes/solanum_lycopersicum_LA2093/rename.genome.fa \
  --outdir    results/TCZZSL20K-LA2093 \
  -profile    standard \
  -resume
