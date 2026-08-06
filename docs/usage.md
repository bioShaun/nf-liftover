# nf-liftover 使用说明

`nf-liftover` 将原来前 4 个手工步骤收敛为一条 Nextflow 命令：输入 `chrom_pos` 格式 ID 文件或 VCF/VCF.GZ/BCF 文件、原始参考基因组 FASTA、目标基因组 FASTA，输出 chain 与 lifted ID/BED 文件或 lifted VCF 文件。

`prepare-probe-dir.sh` 暂不在本流程内，流程结束后仍按原方式手工执行。

## 环境

本服务器上使用 miniforge/mamba：

```bash
MAMBA=/project/software/miniforge3/bin/mamba
$MAMBA run -n nextflow26 nextflow -version
```

当前配置默认使用以下 conda env：

- `nextflow26`：运行 Nextflow 26.x 与 Java。
- `nf-liftover-tools`：运行流程任务，包含 `python`、`pandas`、`typer`、`loguru`、`pyfaidx`、`samtools`、`minimap2`、`transanno`。

在新机器上可从仓库根目录的 `environment.yml` 创建任务环境：

```bash
mamba env create -f environment.yml
```

如果 env 所在目录不同，可用 `--conda_dir` 覆盖；如果 env 名不同，可在 `nextflow.config` 的 `params.conda_envs` 中调整。

查看帮助（无需输入文件）：

```bash
nextflow run . --help
```

## 配置初始化

首次下载或在新的机器上运行时，需要从模板配置文件复制：

```bash
cp nextflow.config.example nextflow.config
```

然后在本地生成的 `nextflow.config` 中根据本机的 conda 路径调整 `conda_dir` 等参数（本地 `nextflow.config` 已加入 `.gitignore` 中，不会提交到仓库）。

## 基本运行

ID 输入模式：

```bash
/project/software/miniforge3/bin/mamba run -n nextflow26 \
  nextflow run /public/scripts/nf-liftover \
  --id        /path/to/TCZZSL20K.id \
  --ref_fa    /path/to/S_lycopersicum_chromosomes.4.00.fa \
  --query_fa  /path/to/rename.genome.fa \
  --outdir    results/TCZZSL20K-LA2093 \
  -profile    standard \
  -resume
```

VCF 输入模式：

```bash
/project/software/miniforge3/bin/mamba run -n nextflow26 \
  nextflow run /public/scripts/nf-liftover \
  --vcf       /path/to/input.vcf.gz \
  --ref_fa    /path/to/source.fa \
  --query_fa  /path/to/target.fa \
  --mapping   /path/to/source-to-target.mapping.tsv \
  --outdir    results/input-liftover \
  -profile    standard \
  -resume
```

HPC 上把 `-profile standard` 改为 `-profile slurm` 或 `-profile slurm_new`。

## 常用参数

- `--id`：输入 ID 文件，每行一个 `chrom_pos` ID；与 `--vcf` 二选一。
- `--vcf`：输入 VCF/VCF.GZ/BCF 文件；与 `--id` 二选一。VCF 模式会保留 FORMAT、样本列、基因型和 INFO 字段。
- `--ref_fa`：原始/参考 FASTA。
- `--query_fa`：目标/新 FASTA。
- `--chain`：可选，已有 chain 文件；提供时跳过 minimap2 与 PAF-to-chain，并校验 chain header 与当前基因组是否匹配。
- `--chain_meta`：可选，先前运行生成的 `chain_meta.yml`；与 `--chain` 联用时额外校验 FASTA sha256。
- `--outdir`：输出目录，默认 `results`。
- `--mapping`：可选，两列 TSV，覆盖自动染色体对应。
- `--pair_strategy`：自动对应策略，`suffix` 或 `order`，默认 `suffix`。
- `--align_mode`：比对模式，`auto`、`whole` 或 `split`，默认 `auto`；`whole` 对应旧 `align_chromosomes-whole.nf`，`split` 对应旧 `align_chromosomes.nf`。
- `--split_threshold`：`align_mode=auto` 时，超过该长度的染色体走 sliding split，默认 `100000000`。
- `--split_size`：split 模式 sliding window 大小，默认 `10000000`。
- `--split_bed`：可选，4 列 `chrom, split_start, split_end, new_chrom` BED；提供后 `<probe>.bed` 与 `<probe>.snpcalling.bed` 输出 split 坐标。
- `--split_genome_fai`：可选，split genome `.fai`，用于排序 split 坐标 BED；未提供时会尝试使用 `split_bed` 同目录下的 `genome.fa.fai`。
- `--flank`：`snpcalling.bed` 两侧扩展长度，默认 `100`。
- `--publish_paf`：是否发布合并后的 `all.paf`，默认 `false`。

## 输出

流程会发布到 `--outdir`：

- `chain/all.chain`（端到端构建或 reuse 校验后发布的副本）
- `chain/chain_meta.yml`（源/目标 FASTA 校验和、染色体映射、校验结果）
- `chain/all.paf`（仅 `--publish_paf true`）
- `software_versions.yml`
- `run_meta.yml`（含 git commit / dirty / manifest 版本）

ID 模式额外输出：

- `liftover/<probe>.id`
- `liftover/<probe>.bed`
- `liftover/<probe>.pos.tsv`
- `liftover/<probe>.snpcalling.bed`

其中 `<probe>` 来自输入 ID 文件名 stem，例如 `TCZZSL20K.id` 对应 `TCZZSL20K.bed`。

VCF 模式额外输出：

- `liftover/<prefix>.vcf.gz`
- `liftover/<prefix>.vcf.gz.tbi`
- `liftover/rejected.<prefix>.vcf.gz`
- `liftover/rejected.<prefix>.vcf.gz.tbi`

其中 `<prefix>` 来自输入 VCF 文件名，去掉 `.vcf.gz`、`.vcf` 或 `.bcf` 后缀。

启用 `--split_bed` 时，`<probe>.id` 与 `<probe>.pos.tsv` 仍记录 liftover 后的目标基因组原始坐标，只有 `<probe>.bed` 与 `<probe>.snpcalling.bed` 转成 split genome 坐标。`--split_bed` 只用于 ID 模式，VCF 模式不做 split-coordinate BED 输出。

## 测试

仓库内包含一个番茄 10 kb smoke fixture。如果 `nextflow26` 环境已安装 nf-test，可跑完整流程回归：

```bash
/project/software/miniforge3/bin/mamba run -n nextflow26 \
  nf-test test tests/tomato_smoke.nf.test
```

当前服务器也可以直接运行 smoke helper,ID 与 VCF 两种模式：

```bash
# ID 模式（默认）
bash tests/data/tomato-smoke/run-smoke.sh id

# VCF liftover 模式
bash tests/data/tomato-smoke/run-smoke.sh vcf

# 复用已有 chain（需先跑 id）
bash tests/data/tomato-smoke/run-smoke.sh chain
```

## 参数校验

流程入口实现了内置 `--help` 与 fail-fast 校验（双层：手写校验先行，随后由 `nf-schema` 插件按 `nextflow_schema.json` 做 schema 级校验）：

- 必填：`--ref_fa`、`--query_fa`，以及 `--id` / `--vcf` 二选一；
- 枚举：`--align_mode`、`--pair_strategy`；
- 范围：`--split_threshold`、`--split_size`、`--flank`、`--max_cpus`；
- 文件存在性：FASTA、FAI、mapping、chain、ID/VCF、split BED；
- 互斥：VCF 模式不可同时提供 `--split_bed` / `--split_genome_fai`。

schema 层（`nf-schema@2.7.2`，在 `nextflow.config` 的 `plugins` 块声明，为启动依赖）额外校验参数类型、`--id` / `--vcf` 的 anyOf 互斥以及 `split_genome_fai → split_bed` 依赖关系。离线环境需预先在 `~/.nextflow/plugins/` 缓存该插件（本机有网时首跑自动拉取），否则 Nextflow 在启动解析阶段即失败。

`nextflow_schema.json` 与 `nextflow.config` / README 中的默认值保持同步，供文档与外部工具使用。
