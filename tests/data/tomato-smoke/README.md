# Tomato Smoke Test Data

Small FASTA fixtures derived from the tomato SL4.0 chromosome 1 sequence.

- `SL4.0ch01.10kb.fa`: first 10 kb from `SL4.0ch01`.
- `LA2093.chr01.10kb.fa`: same 10 kb sequence with the target-style chromosome name `chr01`.
- `chrom_pairs.tsv`: maps `SL4.0ch01` to `chr01`.
- `tomato-smoke.id`: three `chrom_pos` IDs within the 10 kb region.

The target FASTA intentionally uses a homologous renamed copy of the source segment so the full Nextflow smoke test can produce a deterministic chain quickly.
