# nf-liftover 使用说明

`nf-liftover` 将原来前 4 个手工步骤收敛为一条 Nextflow 命令：输入 `chrom_pos` 格式 ID 文件、原始参考基因组 FASTA、目标基因组 FASTA，输出 chain 与 lifted ID/BED 文件。

`prepare-probe-dir.sh` 暂不在本流程内，流程结束后仍按原方式手工执行。

## 环境

建议先进入 Nextflow 环境：

```bash
conda activate nextflow
```

当前配置默认使用以下 conda env：

- `bioinfo`：`seqkit`、`samtools`、`minimap2`、`transanno`
- `py-13`：`python`、`pandas`、`typer`、`loguru`、`pyfaidx`

如果 env 所在目录不同，可用 `--conda_dir` 覆盖；如果 env 名不同，可在 `nextflow.config` 的 `params.conda_envs` 中调整。

## 基本运行

```bash
nextflow run /public/scripts/nf-liftover \
  --id        /path/to/TCZZSL20K.id \
  --ref_fa    /path/to/S_lycopersicum_chromosomes.4.00.fa \
  --query_fa  /path/to/rename.genome.fa \
  --outdir    results/TCZZSL20K-LA2093 \
  -profile    standard \
  -resume
```

HPC 上把 `-profile standard` 改为 `-profile slurm` 或 `-profile slurm_new`。

## 常用参数

- `--id`：输入 ID 文件，每行一个 `chrom_pos` ID。
- `--ref_fa`：原始/参考 FASTA。
- `--query_fa`：目标/新 FASTA。
- `--outdir`：输出目录，默认 `results`。
- `--mapping`：可选，两列 TSV，覆盖自动染色体对应。
- `--pair_strategy`：自动对应策略，`suffix` 或 `order`，默认 `suffix`。
- `--split_threshold`：超过该长度的染色体走 sliding split，默认 `100000000`。
- `--split_size`：大染色体 sliding window 大小，默认 `10000000`。
- `--flank`：`snpcalling.bed` 两侧扩展长度，默认 `100`。

## 输出

流程会发布到 `--outdir`：

- `chain/all.chain`
- `liftover/<probe>.id`
- `liftover/<probe>.bed`
- `liftover/<probe>.pos.tsv`
- `liftover/<probe>.snpcalling.bed`

其中 `<probe>` 来自输入 ID 文件名 stem，例如 `TCZZSL20K.id` 对应 `TCZZSL20K.bed`。

## nf-schema 插件

本流程使用 `nf-schema@2.2.0` 做参数校验。联网环境首次运行会自动下载插件；离线 HPC 可先在有网机器运行一次：

```bash
/project/software/miniforge3/envs/nextflow/bin/nextflow config /public/scripts/nf-liftover
```

然后把 `~/.nextflow/plugins/nf-schema-2.2.0` 同步到 HPC 用户的 `~/.nextflow/plugins/` 下。
