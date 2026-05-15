# Liftover 流程说明（示例）

> Legacy 手工流程参考。新流程请优先使用 [`docs/usage.md`](usage.md) 中的一条 `nextflow run /public/scripts/nf-liftover` 命令。

以 **Solanum lycopersicum (SL4.0)** 坐标的变异位点转换至 **LA2093** 参考基因组为例。

## 流程参数一览

| 维度 | 信息 |
|------|------|
| Input BED/ID | `TCZZSL20K.id` |
| Ref Genome (Q) | `~/genome/solanum_lycopersicum/v4/genome.fa` |
| Target Genome (T) | `~/genome/LA2093/rename.genome.fa` |
| Output Dir | `./TCZZSL20K-LA2093` |
| Main Tools | seqkit, nextflow, transanno |

---

## Step 1：染色体拆分（Genome Splitting）

为支持后续并行或按染色体比对，使用 **seqkit** 按染色体 ID 拆分基因组。

Liftover 工作目录建议：`~/project/standard/genome_alignment/`

### 拆分 Target：LA2093

```bash
conda activate ngs

seqkit split -i \
  --by-id-prefix "la2093." \
  --id-regexp "" \
  -O split-fa/ rename.genome.fa
```

### 拆分 Query：SL4.0

```bash
seqkit split -i \
  --by-id-prefix "" \
  --id-regexp "" \
  -O split-fa/ genome.fa
```

---

## Step 2：准备染色体对应表（Chromosome Mapping）

建立 Query（SL4.0）与 Target（LA2093）的对应列表；**顺序会影响后续比对流程**。

| SL4.0 (Query) | LA2093 (Target) |
|---------------|-----------------|
| SL4.0ch01 | la2093.chr01 |
| SL4.0ch02 | la2093.chr02 |
| SL4.0ch03 | la2093.chr03 |
| SL4.0ch04 | la2093.chr04 |
| SL4.0ch05 | la2093.chr05 |
| SL4.0ch06 | la2093.chr06 |
| SL4.0ch07 | la2093.chr07 |
| SL4.0ch08 | la2093.chr08 |
| SL4.0ch09 | la2093.chr09 |
| SL4.0ch10 | la2093.chr10 |
| SL4.0ch11 | la2093.chr11 |
| SL4.0ch12 | la2093.chr12 |

---

## Step 3：生成 Chain 文件（Alignment Workflow）

使用 **Nextflow** 跑比对流程：内部为 **minimap2 -cx asm5** 全基因组比对，并用 **transanno minimap2chain** 生成 `.chain` 文件。

```bash
source /public/home/zxchen/software/miniconda3/etc/profile.d/conda.sh
conda activate nextflow24
```

### 单条染色体小于 100 Mb

```bash
nextflow run ~/scripts/ngs-utils-v24/workflows/align_chromosomes-whole.nf \
  --chr_pairs sl4-vs-la2093 \
  --species_a_dir v4 \
  --species_b_dir la2093 \
  --ref_fai /Query-genome/genome.fa.fai \
  -profile slurm_new2 \
  -resume
```

### 单条染色体大于等于 100 Mb

```bash
nextflow run ~/scripts/ngs-utils-v24/workflows/align_chromosomes.nf \
  --chr_pairs sl4-vs-la2093 \
  --species_a_dir v4 \
  --species_b_dir la2093 \
  --ref_fai /Query-genome/genome.fa.fai \
  -profile slurm_new2 \
  -resume
```

说明：`--ref_fai` 为 **Query 参考基因组**对应的 `.fai` 路径，请按实际部署修改。

```bash
conda deactivate
cat combined/* > all.paf
```

---

## Step 4：坐标转换（Transanno Liftover）

基于 Step 3 的 `.chain`，将 BED 或 ID 坐标从旧基因组投影到新基因组。

```bash
conda activate py3-13
```

### 不需要拆分

参数含义：待转换 ID 文件、chain、Query FASTA、Target FASTA、输出路径前缀。

```bash
python ~/scripts/tc-pytools/liftover/liftover_by_id.py \
  TCZZSL20K.id \
  results/all.chain \
  ~/genome/solanum_lycopersicum/v4/genome.fa \
  ~/genome/LA2093/rename.genome.fa \
  ./TCZZSL20K-LA2093
```

### 需要拆分（小麦示例）

```bash
python ~/scripts/tc-pytools/liftover/liftover_by_id.py \
  wheat.20K.id \
  results/all.chain \
  ~/genome/nf/iwgsc2.1/genome.fa \
  ~/genome/nf/wheat_ak/genome.fa \
  ./iwgsc2ak \
  --split-bed ~/genome/nf/wheat_ak_split/split.cat.bed \
  --split-genome-fai ~/genome/nf/wheat_ak_split/genome.fa.fai
```

---

## Step 5：结果整理与目录准备

将转换结果整理到统一目录，便于后续分析。

```bash
conda activate py3-13

sh ~/scripts/generic/probe-design/prepare-probe-dir.sh \
  TCZZSL20K-LA2093 \
  solanum_pimpinellifolium_LA2093 \
  TCZZSL20K \
  -n 20
```

### `prepare-probe-dir.sh` 用法摘要

```text
sh ~/scripts/generic/probe-design/prepare-probe-dir.sh <probe_dir> <genome> <probe_id> -n <number>
```

- **`-n`**：拆分 SNP Calling BED 的数量。
