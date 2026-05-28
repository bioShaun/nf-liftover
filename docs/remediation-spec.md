# nf-liftover 修复与优化 Spec

## 背景

当前 `nf-liftover` 已经具备 Nextflow DSL2 typed workflow、`nf-schema` 参数校验、subworkflow 拆分和 smoke fixture。但最近检查发现，仓库里仍混入本机/部署环境配置，部分运行语义会浪费资源或影响结果正确性，后续模块化也缺少明确落地顺序。

本 spec 用于约束下一轮修复：先保证仓库可提交、测试可复现、运行失败行为合理，再逐步迁移到更接近 nf-core 风格的模块化结构。

## 目标

1. 仓库提交 example 配置，不提交个人或单机 `nextflow.config`。
2. 删除手写 help 和 `System.exit(0)`，由 `nf-schema` 管理参数帮助与校验。
3. 只对资源或信号类失败做 retry，避免确定性失败重复消耗资源。
4. 为 retry 后的资源请求设置硬上限。
5. 修复可选输入文件 staging、split BED 跨边界坐标转换等正确性问题。
6. 增加 report、timeline、trace 等运行可观测性输出。
7. 为后续模块化、`ext.args`、per-process versions 和 stub tests 留出结构。

## 非目标

- 不把 `prepare-probe-dir.sh` 纳入本轮流水线。
- 不改变核心用户输入输出契约：输入 ID、ref FASTA、query FASTA，输出 chain、lifted ID/BED、snpcalling BED。
- 不在第一阶段引入 Docker/Singularity 镜像体系；先保留 conda 路径可覆盖机制。
- 不一次性重写所有 subworkflow。模块化应分阶段迁移，避免大范围行为漂移。

## 问题清单

### P0: 配置提交策略错误

`nextflow.config` 当前包含宿主机相关默认值，例如 `conda_dir` 和 env 名映射，却已被 git 跟踪。该文件应作为本地运行配置保留在开发者机器上，不应作为仓库标准配置提交。

要求：

- 提交 `nextflow.config.example`。
- 将 `nextflow.config` 加入 `.gitignore`。
- 从 git index 移除 `nextflow.config`，但不删除本地文件。
- 文档说明首次运行时复制 example：

```bash
cp nextflow.config.example nextflow.config
```

验收：

- `git ls-files nextflow.config` 不再返回文件。
- `git ls-files nextflow.config.example` 返回文件。
- 本地 `nextflow.config` 可继续存在且不显示为待提交变更。

### P0: `System.exit(0)` 和手写 help

`main.nf` 中 `if (params.help)` 分支调用 `System.exit(0)`，会强行终止 JVM，并绕过 Nextflow 正常生命周期。项目已经使用 `nf-schema`，应删除整段手写 help。

要求：

- 删除 `main.nf` 中 `params.help` 条件分支。
- 删除重复维护的 help 文本。
- 依赖 `nf-schema` schema 生成参数帮助和校验。

验收：

- `rg "System.exit|params.help" main.nf` 无结果。
- 参数缺失时仍由 `validateParameters()` 报错。

### P0: `errorStrategy` 过宽

当前配置对所有失败都 retry。语法错误、路径不存在、命令不存在、输入校验失败不应 retry。

要求：

```groovy
process {
    errorStrategy = {
        task.exitStatus in ((130..145) + 104) ? 'retry' : 'terminate'
    }
    maxRetries = 2
}
```

验收：

- `samtools: command not found` 这类 exit 127 不 retry。
- OOM 或信号退出码仍可 retry。

### P0: 缺少资源上限

`large_mem` 使用 `{ 60.GB * task.attempt }`，第三次尝试会请求 180GB，可能超出节点限制，也与 `params.max_memory` 的语义冲突。

要求：

- 使用 `resourceLimits` 或等价 clamp 机制。
- `params.max_cpus`、`params.max_memory`、`params.max_time` 是真实上限，不只是默认值。
- 所有 retry 放大的资源不得超过上限。

建议配置形态：

```groovy
process {
    resourceLimits = [
        cpus: params.max_cpus,
        memory: params.max_memory,
        time: params.max_time
    ]
}
```

如果目标 Nextflow 版本对 `resourceLimits` 行为不满足需求，则引入本地 `check_max()` helper，但必须有测试或 dry-run 验证。

验收：

- `large_mem` 第三次尝试不超过 `params.max_memory`。
- `split_mem` retry 后也不超过配置上限。

### P0: 可选 `.fai` 未作为 `Path` stage

`MAYBE_FAIDX` 中 `fai_hint` 是 `String`。当用户提供 `--ref_fai` 或 `--query_fai` 时，Nextflow 不会 stage 该文件，容器或远程 executor 下可能读取不到。

要求：

- 将可选 `.fai` 输入建模为 `Path?`。
- `main.nf` 中 `params.ref_fai` 和 `params.query_fai` 存在时使用 `file()`。
- process 内优先使用 staged `.fai`；缺失时才运行 `samtools faidx`。

验收：

- 提供 `--ref_fai/--query_fai` 时，`MAYBE_FAIDX` 不重新生成索引。
- 相对路径和绝对路径都能被 stage 到 task work dir。

### P0: split BED 跨边界坐标转换错误

当前 `split_bed_dataframe()` 只按区间 start 所在 split 片段归属，不处理区间跨 split 边界的情况。例如原坐标 `chr01:90-110` 落入 `0-100` 和 `100-200` 两段时，应输出两个 split 区间，而不是单个越界区间。

要求：

- 按 interval overlap 计算 split 区间。
- 对每个 overlap 输出一行。
- 输出坐标 clamp 到 split segment 边界。
- 保持 `.id` 和 `.pos.tsv` 仍为目标基因组原始坐标；只转换 `.bed` 和 `.snpcalling.bed`。

验收：

- 新增单测覆盖跨边界区间。
- 不跨边界的现有 split fixture 输出保持不变。

### P1: `restore_split_paf.py` 参数命名不准

`restore_split_paf.py` 的 `query_fai` 实际传入的是 ref/original FASTA index，用于把 split-window query name 还原到原始 ref chromosome 长度和坐标。

要求：

- 将函数参数、CLI 参数和 help 文本中的 `query_fai` 改为 `ref_fai`。
- 同步调用点和测试名称。

验收：

- `rg "query_fai" bin/restore_split_paf.py subworkflows/align_and_chain.nf tests/test_core_helpers.py` 不再出现旧语义。
- 现有 restore split PAF 单测通过。

### P1: 缺少运行报告输出

应默认生成 Nextflow report、timeline、trace，便于 HPC 排查和性能回顾。

要求：

- 在 example config 中启用：
  - `trace.enabled = true`
  - `report.enabled = true`
  - `timeline.enabled = true`
- 输出路径放在 `${params.outdir}/pipeline_info/`。
- 可选生成 `dag`，但不应阻塞无 graphviz 的环境。

验收：

- smoke run 后 `pipeline_info/` 下存在 trace、report、timeline。

### P1: profiles 配置分散

当前 `profiles {}` 在 `conf/slurm.config`，但作为用户入口配置，profiles 更适合集中在 `nextflow.config.example`。

要求：

- 将 `standard`、`slurm`、`slurm_new` profiles 合并到 `nextflow.config.example`。
- `conf/base.config` 保留资源 tier。
- `conf/conda.config` 保留 env label 映射。
- 如保留 `conf/slurm.config`，只能放 executor 细节，不再定义 profiles。

验收：

- 用户阅读 `nextflow.config.example` 即可看到可用 profiles。
- `nextflow config -profile standard` 能解析完整配置。

### P2: 缺少 `modules.config` 和 `ext.args`

minimap2、transanno 等 CLI 参数硬编码在 process script 中，不便于按 profile 或单个 process 调参。

要求：

- 新增 `conf/modules.config`。
- 在 process 中使用 `task.ext.args ?: ''` 和必要的 `task.ext.prefix`。
- 先覆盖以下进程：
  - `BUILD_QUERY_MMI`
  - `ALIGN_WHOLE_CHROMOSOME`
  - `ALIGN_SPLIT_WINDOW`
  - `PAF_TO_CHAIN`
  - `LIFTOVER_BY_ID`

验收：

- 修改 minimap2 参数只需要改 config，不需要改 `.nf`。
- 默认行为与现有命令等价。

### P2: process 仍集中在 subworkflow 文件中

当前 `subworkflows/*.nf` 同时定义 process 和 workflow，不利于复用和测试。

要求：

- 新建 `modules/local/`。
- 每个 process 一个 module 文件。
- `subworkflows/` 只负责 include modules、组合 channels、emit outputs。
- 分批迁移，迁移一组跑一次测试。

建议顺序：

1. `MAYBE_FAIDX`、`DERIVE_CHROM_PAIRS`
2. `BUILD_ALIGN_PAIRS`、`EXTRACT_CHROM_FASTAS`
3. minimap2 和 PAF/chain 相关进程
4. `LIFTOVER_BY_ID`
5. `COLLECT_SOFTWARE_VERSIONS` 或版本合并进程

验收：

- 每次迁移后 `nf-test` smoke 通过。
- subworkflow 文件中不再出现 `process` 定义。

### P2: 版本收集不标准

当前 `COLLECT_SOFTWARE_VERSIONS` 要求一个进程环境同时包含 python、samtools、seqkit、minimap2、transanno。该模式会导致胖环境，并和模块化原则冲突。

要求：

- 每个分析 process 输出自己的 `versions.yml`。
- workflow 末尾合并所有 `versions.yml`。
- 合并进程只依赖轻量 Python 或 shell 环境。

验收：

- `software_versions.yml` 包含每个实际执行工具版本。
- 不再需要 `tool_py_ngs` 同时装所有工具来查询版本。

### P2: 缺少 `stub:`

没有 stub 时，拓扑测试必须依赖真实 `samtools/minimap2/transanno`。这会让本地和 CI 的 smoke test 难以稳定。

要求：

- 为核心 process 增加 `stub:`。
- stub 输出满足下游 channel shape 和文件名契约。
- nf-test 增加一组 `-stub-run` 拓扑测试。

验收：

- `nf-test` 可以在没有真实生信工具的环境跑拓扑测试。
- 真实 smoke test 仍保留，用于完整回归。

## 实施阶段

### 阶段 1: 配置与测试基线

交付：

- `nextflow.config.example`
- `.gitignore` 忽略本地 `nextflow.config`
- `nf-test.config` 使用 example 或专用 test config
- `docs/usage.md` 更新本地配置说明

验证：

```bash
git ls-files nextflow.config
git ls-files nextflow.config.example
mamba run -n probe-design pytest tests/test_core_helpers.py
```

### 阶段 2: P0 运行语义修复

交付：

- 删除 `System.exit(0)` help 分支
- 限定 `errorStrategy`
- 增加资源上限
- 修复 `.fai` Path staging
- 修复 split BED 跨边界转换

验证：

```bash
mamba run -n probe-design pytest tests/test_core_helpers.py
mamba run -n nextflow nf-test test tests/tomato_smoke.nf.test
```

### 阶段 3: P1 可维护性修复

交付：

- `restore_split_paf.py` 参数改名
- report/timeline/trace 输出
- profiles 合并进 example config

验证：

```bash
mamba run -n nextflow nextflow config -profile standard
mamba run -n nextflow nf-test test tests/tomato_smoke.nf.test
```

### 阶段 4: P2 模块化和快速测试

交付：

- `conf/modules.config`
- `modules/local/`
- per-process `versions.yml`
- `stub:` blocks
- stub topology nf-test

验证：

```bash
mamba run -n nextflow nf-test test tests/tomato_smoke.nf.test
mamba run -n nextflow nf-test test tests/tomato_smoke.nf.test --profile standard
```

如果 nf-test 对 `-stub-run` 的参数传递方式需要调整，以实际 nf-test 版本命令为准。

## 最终验收标准

- `nextflow.config` 不再被 git 跟踪。
- `nextflow.config.example` 可作为新机器初始化模板。
- Python 单元测试全过。
- Nextflow smoke test 在指定测试环境全过。
- 缺失命令、语法错误、输入校验失败不再 retry。
- OOM/信号类退出码仍会 retry，但资源请求不超过上限。
- 提供 `.fai` 时可被 Nextflow stage 并复用。
- split BED 跨边界区间输出正确拆分。
- 运行结果包含 trace、report、timeline。
- 后续 process 模块化可以小步迁移，不改变用户命令和输出契约。
