# nf-liftover

`nf-liftover` 是一个基于 Nextflow DSL2 的端到端基因组坐标转换流程。它从源基因组和目标基因组生成染色体级 PAF/chain 文件，并使用 chain 将位点 ID 或 VCF 变异投影到目标基因组。

流程支持：

- `chrom_pos` 格式的位点 ID 文件和 VCF/VCF.GZ/BCF 两种输入模式；
- 按染色体名称后缀、染色体顺序或显式映射表建立染色体对应关系；
- 整条染色体比对、滑窗拆分比对，以及按染色体长度自动选择比对模式；
- 本地和 Slurm 执行；
- 自动生成 FASTA 索引、运行报告、软件版本记录和可恢复的 Nextflow 工作目录。

## 工作流程

```text
源 FASTA ─┐
          ├─ FASTA 索引 ─ 染色体配对 ─ minimap2 ─ PAF ─ chain ─┬─ ID liftover
目标 FASTA ┘                                                    └─ VCF liftover
```

当 `--align_mode auto` 时，长度小于 `--split_threshold` 的源染色体采用整条比对，其余染色体按 `--split_size` 拆分后并行比对，再恢复为原始坐标并合并 PAF。

## 运行要求

- Linux 和 Bash；
- Nextflow `>=26.0.0` 及兼容的 Java；
- Conda 或 Mamba；
- 任务环境中的 `samtools`、`seqkit`、`minimap2`、`transanno`、Python、pandas、Typer、Loguru 和 pyfaidx。

项目默认从 `/project/software/miniforge3` 查找 Conda，并使用 `nf-liftover-tools` 环境。若本机路径或环境名不同，请修改本地配置。

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

## 输入

以下参数始终必填：

- `--ref_fa`：源坐标所在的参考基因组 FASTA；
- `--query_fa`：需要投影到的目标基因组 FASTA；
- `--id` 或 `--vcf`：二选一，不能同时提供。

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

在集群上将 `standard` 替换为 `slurm`；`slurm_new` 使用 `cae` 队列。队列和资源上限可分别在 `conf/slurm.config` 与本地 `nextflow.config` 中调整。

## 主要参数

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--ref_fai` | 自动生成 | 可选的源 FASTA `.fai` |
| `--query_fai` | 自动生成 | 可选的目标 FASTA `.fai` |
| `--mapping` | 自动配对 | 两列染色体映射 TSV |
| `--pair_strategy` | `suffix` | 自动配对策略：`suffix` 或 `order` |
| `--align_mode` | `auto` | `auto`、`whole` 或 `split` |
| `--split_threshold` | `100000000` | `auto` 模式切换到滑窗比对的染色体长度 |
| `--split_size` | `10000000` | 滑窗大小 |
| `--split_bed` | 无 | ID 模式下将 BED 输出转换为 split-genome 坐标的四列 BED |
| `--split_genome_fai` | 自动推断 | split genome 的 `.fai`；默认查找 `split_bed` 同目录的 `genome.fa.fai` |
| `--flank` | `100` | `snpcalling.bed` 向两侧扩展的碱基数 |
| `--outdir` | `results` | 发布结果目录 |
| `--max_cpus` | `16` | 单任务 CPU 上限 |
| `--max_memory` | `60.GB` | 单任务内存上限 |
| `--max_time` | `48.h` | 单任务运行时间上限 |

完整参数定义见 [`nextflow_schema.json`](nextflow_schema.json)。

## 输出

所有模式都会生成：

```text
<outdir>/
├── chain/
│   ├── all.paf
│   └── all.chain
├── pipeline_info/
│   ├── execution_report.html
│   ├── execution_timeline.html
│   └── execution_trace.txt
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
```

## 项目结构

```text
main.nf                 主工作流入口
subworkflows/           基因组准备、比对和 liftover 子工作流
modules/local/          Nextflow process 模块
bin/                    坐标、PAF 和结果处理脚本
conf/                   本地、测试和 Slurm 配置
tests/                  Python 单元测试与 nf-test 测试
examples/               可运行示例及染色体映射
docs/                   使用说明、历史流程和改进设计文档
```

更详细的参数和 split-genome 说明见 [`docs/usage.md`](docs/usage.md)，番茄示例见 [`examples/sl4-vs-la2093/`](examples/sl4-vs-la2093/)。`docs/pipe.md` 记录的是旧版手工流程，仅供迁移和结果核对使用。

## 注意事项

- 当前 `conf/conda.config` 的 `beforeScript` 会将 Miniforge base 的 `bin` 放在任务环境之前。在从 `nextflow26` 环境启动流程时，Python 任务可能错误使用 base Python。正式运行前应先修正该 PATH 优先级；详见下方“已知优化项”。
- `ref` 表示位点当前所在的源基因组，`query` 表示目标基因组；运行前应确认方向，避免生成反向 chain。
- `suffix` 自动配对依赖染色体名中可识别且唯一的数字后缀；命名不规则时应显式提供 `--mapping`。
- VCF 模式保留 FORMAT、样本、基因型和 INFO 字段，无法转换的记录写入 `rejected.*.vcf.gz`。
- `-resume` 依赖原工作目录，清理 `work/` 后无法复用历史任务。
- `prepare-probe-dir.sh` 不属于本项目，若下游流程需要，仍需在 liftover 完成后单独执行。

## 已知优化项

当前最优先的改进是修正 `conf/conda.config` 的任务环境激活，确保 `python`、`samtools`、`minimap2` 和 `transanno` 都来自声明的 `nf-liftover-tools` 环境。随后建议补充可复现的环境定义或容器、真正的 `--help` 和参数校验、CI 中的 ID/VCF smoke test，以及允许复用已有 chain 文件的独立 liftover 入口。完整证据、优先级和验收标准见 [`docs/optimization-assessment.md`](docs/optimization-assessment.md)。
