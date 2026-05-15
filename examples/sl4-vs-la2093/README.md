# SL4.0 到 LA2093 示例

这个示例使用已就位的数据：

- `/public/data/genomes/solanum_lycopersicum_LA2093/TCZZSL20K.id`
- `/public/data/genomes/solanum_lycopersicum_LA2093/S_lycopersicum_chromosomes.4.00.fa`
- `/public/data/genomes/solanum_lycopersicum_LA2093/rename.genome.fa`

运行：

```bash
conda activate nextflow
bash examples/sl4-vs-la2093/run.sh
```

输出目录为 `results/TCZZSL20K-LA2093`。流程结束后，如需继续整理 probe 目录，仍按现有方式手工运行 `prepare-probe-dir.sh`。

`chrom_pairs.tsv` 显式映射 `SL4.0ch01-12` 到 `chr01-12`。SL4.0 的 `SL4.0ch00` 在 LA2093 目标 FASTA 中没有对应染色体，因此真实示例不使用自动 suffix 推断。
