# SL4.0 到 LA2093 示例

这个示例使用已就位的数据：

- `/public/data/genomes/solanum_lycopersicum_LA2093/TCZZSL20K.id`
- `/public/data/genomes/solanum_lycopersicum_LA2093/S_lycopersicum_chromosomes.4.00.fa`
- `/public/data/genomes/solanum_lycopersicum_LA2093/rename.genome.fa`

> ⚠️ **注意（2026-08-05 审计，见 [drift-audit-2026-08-05.md](../../docs/drift-audit-2026-08-05.md)）**：上述 `/public/data/genomes/solanum_lycopersicum_LA2093/` 目录在当前主机上不存在（复数和单数 genome 路径均已核实缺失）。运行前需自备 `--id` / `--ref_fa` / `--query_fa` 指向实际可用的文件，或恢复 LA2093 数据集。`run.sh` 中的路径仅作格式示例。

运行：

```bash
conda activate nextflow26
bash examples/sl4-vs-la2093/run.sh
```

输出目录为 `results/TCZZSL20K-LA2093`。流程结束后，如需继续整理 probe 目录，仍按现有方式手工运行 `prepare-probe-dir.sh`。

`chrom_pairs.tsv` 显式映射 `SL4.0ch01-12` 到 `chr01-12`。SL4.0 的 `SL4.0ch00` 在 LA2093 目标 FASTA 中没有对应染色体，因此真实示例不使用自动 suffix 推断。
