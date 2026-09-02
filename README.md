# nf-liftover

`nf-liftover` 是一个基于 Nextflow DSL2 的端到端基因组坐标转换流程。它从源基因组和目标基因组生成染色体级 PAF/chain 文件，并使用 chain 将位点 ID 或 VCF 变异投影到目标基因组。

流程支持：

- `chrom_pos` 格式的位点 ID 文件和 VCF/VCF.GZ/BCF 两种输入模式；
- 按染色体名称后缀、染色体顺序或显式映射表建立染色体对应关系；
- 整条染色体比对、滑窗拆分比对，以及按染色体长度自动选择比对模式；
- 复用已有 `--chain`，跳过昂贵的基因组比对；
- 本地和 Slurm 执行；
- 自动生成 FASTA 索引、运行报告、软件版本与可追溯的运行元数据。

## 工作流程

```text
源 FASTA ─┐
          ├─ FASTA 索引 ─ 染色体配对 ─ 比对(minimap2/mm2plus) ─ PAF ─ chain ─┬─ ID liftover
目标 FASTA ┘          （或 --chain 跳过比对）                                └─ VCF liftover
```

当 `--align_mode auto` 时，长度小于 `--split_threshold` 的源染色体采用整条比对，其余染色体按 `--split_size` 拆分后并行比对，再恢复为原始坐标并合并 PAF。

## 运行要求

- Linux 和 Bash；
- Nextflow `>=26.0.0` 及兼容的 Java；
- Conda 或 Mamba；
- 任务环境中的 `samtools`、`minimap2`、`transanno`、Python、pandas、Typer、Loguru 和 pyfaidx；
- 可选：`mm2plus`（经 `--aligner mm2plus` 启用），需额外 conda 环境，环境名由 `params.aligner_envs` 指定。

### 可复现环境

仓库提供锁定用途的 [`environment.yml`](environment.yml)。在新机器上：

```bash
mamba env create -f environment.yml
# 或更新已有环境
mamba env update -f environment.yml --prune
```

Nextflow 启动环境可单独安装（本服务器常用 `nextflow26`）。任务默认使用 `/project/software/miniforge3/envs/nf-liftover-tools`；若路径不同，请修改本地配置中的 `params.conda_dir` / `params.conda_envs`。

## 配置

复制配置模板：

```bash
cp nextflow.config.example nextflow.config
```

然后按实际环境修改 `nextflow.config` 中的 `params.conda_dir` 和 `params.conda_envs`。`nextflow.config` 是本地配置文件，已被 Git 忽略；可提交的默认配置维护在 `nextflow.config.example`。

检查 Nextflow：

```bash
/project/software/miniforge3/bin/mamba run -n nextflow26 nextflow -version
```

查看参数帮助（无需输入文件）：

```bash
nextflow run . --help
```

## 输入

以下参数始终必填：

- `--ref_fa`：源坐标所在的参考基因组 FASTA；
- `--query_fa`：需要投影到的目标基因组 FASTA；
- `--id` 或 `--vcf`：二选一，不能同时提供。

可选 `--chain`：提供已有 chain 时跳过 minimap2 / PAF-to-chain，只做基因组准备与 liftover。

ID 文件每行一个位点，使用最后一个下划线分隔染色体和 1-based 位置。例如：

```text
SL4.0ch01_1024
SL4.0ch01_2048
```

显式染色体映射表为无表头的两列 TSV，第一列是源染色体，第二列是目标染色体：

```text
SL4.0ch01	la2093.chr01
SL4.0ch02	la2093.chr02
```

未提供 `--mapping` 时，流程使用 `--pair_strategy suffix`（默认）或 `order` 自动配对。

## 快速开始

### ID 模式

```bash
/project/software/miniforge3/bin/mamba run -n nextflow26 \
  nextflow run . \
  --id /path/to/panel.id \
  --ref_fa /path/to/source.fa \
  --query_fa /path/to/target.fa \
  --outdir results/panel-liftover \
  -profile standard \
  -resume
```

### VCF 模式

```bash
/project/software/miniforge3/bin/mamba run -n nextflow26 \
  nextflow run . \
  --vcf /path/to/input.vcf.gz \
  --ref_fa /path/to/source.fa \
  --query_fa /path/to/target.fa \
  --mapping /path/to/source-to-target.tsv \
  --outdir results/vcf-liftover \
  -profile standard \
  -resume
```

### 复用已有 chain

```bash
/project/software/miniforge3/bin/mamba run -n nextflow26 \
  nextflow run . \
  --id /path/to/panel.id \
  --ref_fa /path/to/source.fa \
  --query_fa /path/to/target.fa \
  --chain /path/to/all.chain \
  --chain_meta /path/to/chain_meta.yml \
  --outdir results/panel-liftover-reuse \
  -profile standard \
  -resume
```

复用模式下流程会：

1. 校验 chain header 中的源/目标染色体名称与长度是否匹配当前 `--ref_fa` / `--query_fa`；
2. 检测明显反向的 chain（tName 落在 query、qName 落在 ref）；
3. 若提供 `--chain_meta`，再核对 sidecar 中的 FASTA sha256；
4. 将校验通过的 chain **复制发布**到 `<outdir>/chain/all.chain`，并写出新的 `chain_meta.yml`。

在集群上将 `standard` 替换为 `slurm`；`slurm_new` 使用 `cae` 队列。队列和资源上限可分别在 `conf/slurm.config` 与本地 `nextflow.config` 中调整。

## 主要参数

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--ref_fai` | 自动生成 | 可选的源 FASTA `.fai` |
| `--query_fai` | 自动生成 | 可选的目标 FASTA `.fai` |
| `--mapping` | 自动配对 | 两列染色体映射 TSV |
| `--chain` | 无 | 已有 chain；提供时跳过比对（校验 header 染色体名/长度/方向） |
| `--chain_meta` | 无 | 可选 sidecar；与 `--chain` 联用时校验源/目标 FASTA sha256 |
| `--pair_strategy` | `suffix` | 自动配对策略：`suffix` 或 `order` |
| `--aligner` | `minimap2` | 比对器：`minimap2`（默认）或 `mm2plus` |
| `--align_mode` | `auto` | `auto`、`whole` 或 `split` |
| `--split_threshold` | `100000000` | `auto` 模式切换到滑窗比对的染色体长度 |
| `--split_size` | `10000000` | 滑窗大小 |
| `--split_bed` | 无 | ID 模式下将 BED 输出转换为 split-genome 坐标的四列 BED |
| `--split_genome_fai` | 自动推断 | split genome 的 `.fai`；默认查找 `split_bed` 同目录的 `genome.fa.fai` |
| `--flank` | `100` | `snpcalling.bed` 向两侧扩展的碱基数 |
| `--publish_paf` | `false` | 是否将合并后的 `all.paf` 发布到结果目录 |
| `--outdir` | `results` | 发布结果目录 |
| `--max_cpus` | `16` | 单任务 CPU 上限 |
| `--max_memory` | `60.GB` | 单任务内存上限（含 retry 放大后的硬上限） |
| `--max_time` | `48.h` | 单任务运行时间上限 |

窗口数远多于可用并发槽位时选 `minimap2`；窗口数少于可用槽位（少染色体 liftover、独占大节点）时选 `mm2plus`。选型依据见 [`docs/wheat-split-benchmark-spec.md`](docs/wheat-split-benchmark-spec.md) 的最终结论表。

完整参数定义见 [`nextflow_schema.json`](nextflow_schema.json)。入口在提交 process 前会校验必填项、枚举、数值范围、文件存在性，以及 ID/VCF 与 split 参数互斥关系。

## 输出

所有模式都会生成：

```text
<outdir>/
├── chain/
│   ├── all.chain           # 新建或复用后发布的 chain 副本
│   ├── chain_meta.yml      # 源/目标 FASTA 校验和、染色体映射、校验结果
│   └── all.paf             # 仅当 --publish_paf true（对齐构建模式）
├── pipeline_info/
│   ├── execution_report.html
│   ├── execution_timeline.html
│   └── execution_trace.txt
├── run_meta.yml            # 输入、参数、git commit/dirty、manifest 版本
└── software_versions.yml
```

ID 模式在 `<outdir>/liftover/` 下额外生成：

- `<panel>.id`：转换后的 `chrom_pos` ID；
- `<panel>.bed`：转换后的 BED；
- `<panel>.pos.tsv`：目标染色体、位置和原始 ID；
- `<panel>.snpcalling.bed`：扩展并合并后的 SNP calling 区间。

VCF 模式在 `<outdir>/liftover/` 下额外生成：

- `<prefix>.vcf.gz` 及其 `.tbi` 索引；
- `rejected.<prefix>.vcf.gz` 及其 `.tbi` 索引。

启用 `--split_bed` 时，仅 ID 模式的 `.bed` 和 `.snpcalling.bed` 转换为 split-genome 坐标；`.id` 和 `.pos.tsv` 仍保留目标基因组原始坐标。

## 测试

运行 Python 单元测试：

```bash
/project/software/miniforge3/bin/mamba run -n nf-liftover-tools \
  python -m unittest discover -s tests -p 'test_*.py' -v
```

安装并确保 `nf-test` 在 `PATH` 中后，运行端到端 smoke test：

```bash
nf-test test tests/tomato_smoke.nf.test
```

也可以直接运行仓库中的番茄 10 kb fixture：

```bash
bash tests/data/tomato-smoke/run-smoke.sh id
bash tests/data/tomato-smoke/run-smoke.sh vcf
bash tests/data/tomato-smoke/run-smoke.sh chain   # 需先跑 id 生成 chain
```

CI 定义见 [`.github/workflows/ci.yml`](.github/workflows/ci.yml)：单元测试、配置解析、`--help`、参数校验失败，以及 ID/VCF（含 chain 复用）smoke。

## 项目结构

```text
main.nf                 主工作流入口
subworkflows/           基因组准备、比对和 liftover 子工作流
modules/local/          Nextflow process 模块
bin/                    坐标、PAF 和结果处理脚本
conf/                   本地、测试和 Slurm 配置
environment.yml         可复现任务 Conda 环境
tests/                  Python 单元测试与 nf-test 测试
examples/               可运行示例及染色体映射
docs/                   使用说明、历史流程和改进设计文档
```

更详细的参数和 split-genome 说明见 [`docs/usage.md`](docs/usage.md)，番茄示例见 [`examples/sl4-vs-la2093/`](examples/sl4-vs-la2093/)。`docs/pipe.md` 记录的是旧版手工流程，仅供迁移和结果核对使用。优化优先级与验收标准见 [`docs/optimization-assessment.md`](docs/optimization-assessment.md)。

## 注意事项

- `ref` 表示位点当前所在的源基因组，`query` 表示目标基因组；运行前应确认方向，避免生成反向 chain。复用 `--chain` 时请对照 `chain_meta.yml` 中的 FASTA 校验和与染色体映射。
- `suffix` 自动配对依赖染色体名中可识别且唯一的数字后缀；命名不规则时应显式提供 `--mapping`。
- VCF 模式保留 FORMAT、样本、基因型和 INFO 字段，无法转换的记录写入 `rejected.*.vcf.gz`。
- `-resume` 依赖原工作目录，清理 `work/` 后无法复用历史任务。
- 资源请求受 `max_cpus` / `max_memory` / `max_time` 硬限制；`large_mem` 从 20 GB 起随 retry 放大（20→40→60），不会超过 `max_memory`。
- `prepare-probe-dir.sh` 不属于本项目，若下游流程需要，仍需在 liftover 完成后单独执行。
