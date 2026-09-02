# nf-liftover 优化评估

评估日期：2026-07-15

## 目的与范围

本文基于当前仓库的 Nextflow 配置、工作流模块、Python 辅助脚本、测试和一次本地 smoke 运行进行评估。目标是记录当前可验证的问题、优化优先级、推荐方案和验收标准，作为后续 Issue 与实施计划的依据。

本文不替代 [`upgrade-plan.md`](upgrade-plan.md)。后者记录历史升级目标，本文描述当前实现仍需处理的事项。

## 当前状态

- Nextflow 配置可由 Nextflow 26.04.4 正常解析。
- Python 单元测试共 21 项，使用 `nf-liftover-tools` 环境运行时全部通过。
- 缺少 `--ref_fa` 等必填参数时，入口会以非零状态退出并给出错误信息。
- ID 模式 smoke 流程可以完成 FASTA 索引、染色体配对、minimap2 比对和 chain 生成，但在 `LIFTOVER_BY_ID` 阶段因任务使用了错误的 Python 环境而失败。
- 仓库包含 ID、VCF、whole、split 和 split-coordinate 等 nf-test 场景，`nextflow26` 环境现已安装 `nf-test 0.9.5`（2026-08-05 审计核实，见 [drift-audit-2026-08-05.md](drift-audit-2026-08-05.md)）。

## 优先级总览

| 优先级 | 优化项 | 类型 | 主要收益 |
| --- | --- | --- | --- |
| P0 | 修复任务环境 PATH 优先级 | 正确性 | 恢复端到端可运行性 |
| P1 | 提供可复现的依赖环境 | 可复现性 | 避免共享环境漂移导致回归 |
| P1 | 建立自动化端到端测试 | 质量保障 | 在合并前发现环境和流程问题 |
| P1 | 完善帮助与参数边界校验 | 用户体验/正确性 | 更早、更明确地阻止错误输入 |
| P1 | 支持复用已有 chain | 性能/架构 | 多批次 liftover 避免重复比对 |
| P2 | 优化大型基因组的存储与资源 | 性能/成本 | 降低 PAF 存储和集群资源消耗 |
| P2 | 改进滑窗边界策略 | 结果质量 | 减少窗口边界附近的潜在漏比对 |
| P2 | 提升报告覆盖与运行可观测性 | 运维 | 避免重复运行时报告写入失败 |
| P3 | 拆分 Python 辅助脚本职责 | 可维护性 | 降低修改和测试成本 |

## 详细评估

### P0：修复任务环境 PATH 优先级

#### 证据

`conf/conda.config` 中三个工具 label 的 `beforeScript` 都执行以下逻辑：

```groovy
export PATH=${params.conda_dir}/condabin:${params.conda_dir}/bin:$PATH
```

Nextflow 生成的任务启动脚本先执行 `beforeScript`，再激活 `nf-liftover-tools`。从 `nextflow26` 环境启动时，手工加入的 Miniforge base 路径不会被 Conda 激活逻辑移除，因此仍排在工具环境之前。

运行时观测结果：

```text
CONDA_PREFIX=/project/software/miniforge3/envs/nf-liftover-tools
command -v python=/project/software/miniforge3/bin/python
```

任务看似激活了正确环境，实际执行的却是 base Python 3.13。该 Python 加载 base 环境中的旧版 pyfaidx，并因缺少 `pkg_resources` 失败。直接执行 `nf-liftover-tools/bin/python` 时，pyfaidx 可以正常导入且全部单元测试通过，说明工具环境本身不是根因。

#### 建议

- 删除会把 base `bin` 放到 PATH 首位的逻辑，只保留 Nextflow 所需的 Conda 初始化方式；或通过启动环境设置 `NXF_CONDA_CACHEDIR`，让任务完全由 Nextflow 管理环境。
- 在每类工具任务中验证 `python`、`samtools`、`minimap2` 和 `transanno` 的实际路径与版本来自声明环境。
- 不要通过修补 base 环境的 `pkg_resources` 掩盖问题，因为其他工具仍可能从错误环境解析。

#### 验收标准

- ID 和 VCF smoke test 均成功完成。
- 任务日志中 `python`、`samtools`、`minimap2` 和 `transanno` 的路径均来自 `nf-liftover-tools`。
- 从 base shell、`conda activate nextflow26` 和 `conda run -n nextflow26` 三种入口启动时行为一致。

### P1：提供可复现的依赖环境

#### 现状与影响

流程依赖固定路径下的共享 Conda 环境，但仓库没有 `environment.yml`、显式 lock 文件或容器定义。共享环境升级后，即使代码未变化，Python 与生物信息工具版本也可能改变。此次 base pyfaidx/setuptools 不兼容就是环境漂移影响运行的实例。

#### 建议

- 提供版本固定的 Conda 环境文件，并在 CI 中从该文件创建环境。
- 对生产/HPC 使用场景，优先提供 Apptainer/Singularity 镜像；保留 Conda profile 作为本地开发入口。
- 在 `software_versions.yml` 中继续记录实际版本，并将关键版本纳入 smoke test 结果。

#### 验收标准

- 在没有预创建 `nf-liftover-tools` 的新环境中，可以按仓库文件完成安装并运行 smoke test。
- 同一版本的环境定义可重复生成一致的关键工具版本。

### P1：建立自动化端到端测试

#### 现状与影响

现有测试覆盖面较好，但没有仓库级 CI 配置，`nf-test` 也不在文档默认的 `nextflow26` 环境内。Python 单元测试无法发现 Nextflow 任务激活了错误解释器，因此 P0 问题只能在实际流程运行时暴露。

#### 建议

- CI 至少运行 Python 单元测试、Nextflow 配置解析和番茄 fixture 的 ID/VCF smoke test。
- whole、split、自动染色体配对和 split-coordinate 场景可以作为较慢的定期或合并前测试。
- 将 `nf-test` 版本纳入可复现环境，不依赖机器上偶然存在的可执行文件。

#### 验收标准

- 每次合并请求都能自动验证 ID 与 VCF 端到端输出。
- CI 对任务环境路径错误、输出缺失和非零 process 状态均能失败。

### P1：完善帮助与参数边界校验

#### 现状与影响

`main.nf` 当前只校验 `ref_fa`、`query_fa` 以及 `id`/`vcf` 二选一。`nextflow_schema.json` 定义了 `help`、枚举和数值范围，但入口没有实现 `--help`，也没有统一使用 schema 校验。

运行时仍可能较晚才发现以下问题：

- FASTA、FAI、mapping 或 split BED 路径不存在；
- `align_mode`、`pair_strategy` 值无效；
- `split_size`、`split_threshold` 或资源上限不合理；
- VCF 模式同时提供只适用于 ID 模式的 split-coordinate 参数。

#### 建议

- 实现无需必填输入即可执行的 `--help`。
- 统一 schema 与运行时校验来源，避免默认值和枚举在多个文件中漂移。
- 对文件存在性、数值范围和互斥参数进行 fail-fast 校验。

#### 验收标准

- `nextflow run . --help` 退出状态为 0，且不要求输入文件。
- 所有无效枚举、范围、缺失文件和冲突参数在提交 process 前失败。
- schema、配置和 README 中的默认值保持一致。

### P1：支持复用已有 chain

#### 现状与影响

顶层工作流总是先运行 `PREPARE_GENOMES` 和 `ALIGN_AND_CHAIN`，再进入 ID 或 VCF liftover。对相同基因组组合的多个 panel/VCF 批次，用户只能依赖 Nextflow cache；更换工作目录、清理 `work/` 或调整无关参数后可能重新执行昂贵的基因组比对。

#### 建议

- 增加可选 `--chain` 参数。提供时跳过 alignment，只校验 chain 与 FASTA 方向后执行 liftover。
- 长期可将 “build-chain” 和 “liftover” 作为两个显式入口或子工作流，同时保留当前端到端默认入口。
- 为 chain 输出记录源/目标 FASTA 校验和、染色体映射和 minimap2 参数，防止复用方向或版本错误的 chain。

#### 验收标准

- 使用已有 chain 时不提交任何 minimap2 或 PAF-to-chain 任务。
- 相同输入通过端到端模式和 `--chain` 模式得到一致的 liftover 输出。
- chain 元数据能明确标识源基因组、目标基因组和染色体映射。

### P2：优化大型基因组的存储与资源

#### 现状与影响

- `all.paf` 总是复制到结果目录，大型基因组可能产生明显的重复存储和发布开销。
- whole 模式任务统一使用 `large_mem`，CPU 和内存按全局上限申请，没有根据染色体长度调整。
- `large_mem` 的内存表达式随重试次数增长，但同时受 `max_memory` 限制，默认配置下重试无法获得更多内存。
- 所有染色体 PAF 通过 `toSortedList` 汇总后才进入合并步骤，形成全局等待点。

#### 建议

- 增加 `--publish_paf` 开关，生产默认只发布 chain，调试时保留 PAF。
- 按染色体长度或历史 trace 动态设置 whole alignment 资源。
- 明确重试的资源增长策略与 `resourceLimits` 的关系。
- 对极大规模染色体集合评估流式或分层合并 PAF，减少全局聚合压力。

#### 验收标准

- 关闭 PAF 发布时，结果目录不包含 `all.paf`，但 chain 与 liftover 结果不变。
- 资源申请不会超过配置上限，重试时的资源变化符合文档说明。
- 使用代表性大型基因组记录基线运行时间、峰值内存和结果目录大小，并验证优化后没有结果差异。

### P2：改进滑窗边界策略

#### 现状与影响

split 模式使用不重叠窗口：下一窗口从前一窗口末端之后开始。跨窗口边界的较长变异、重复区域或复杂比对可能缺少足够上下文。当前 smoke fixture 验证坐标恢复，但不足以评估边界附近的生物学结果完整性。

#### 建议

- 引入可配置的 `split_overlap`，在恢复 PAF 后按明确规则去重重叠命中。
- 增加跨窗口边界的合成测试数据，并与 whole 模式输出比较。
- 在引入 overlap 前先用代表性大型基因组评估漏比对比例，避免无证据增加计算量。

#### 验收标准

- 边界测试中的预期比对不会因窗口切分丢失。
- 合并后的 PAF/chain 不包含由 overlap 引入的重复或坐标冲突。
- overlap 的运行时间与结果质量影响有可复现的基准记录。

### P2：提升报告覆盖与运行可观测性

#### 现状与影响

trace 设置了 `overwrite = true`，但 report 和 timeline 没有相同设置。历史日志已经出现因目标 HTML 存在而无法渲染报告的警告，常见于对同一 `outdir` 使用 `-resume`。

#### 建议

- 明确报告覆盖策略：允许覆盖，或按运行 ID/时间戳保存。
- 在流程结束摘要中显示成功、失败、拒绝记录数量，以及 chain/PAF 的关键统计。
- 将关键输入、参数和 Git revision 写入机器可读的运行元数据。

#### 验收标准

- 对同一 `outdir` 连续执行和 `-resume` 时不会丢失 report/timeline。
- 每次运行都能追溯输入、参数、代码版本和软件版本。

### P3：拆分 Python 辅助脚本职责

#### 现状与影响

`bin/liftover_by_id.py` 同时承担 CLI、输入解析、transanno 调用、DataFrame 处理、split-coordinate 转换和多种文件输出。当前测试覆盖较好，但继续增加 rescue、统计或新输出格式会提高修改风险。

#### 建议

- 保持 CLI 入口轻量，将坐标转换、split 映射和输出写入拆成边界清晰的模块。
- 先锁定当前输出行为，再进行机械拆分；避免在重构时同时改变坐标语义。
- 统一 `liftover_by_id.py` 与 `merge_blast_rescue.py` 的公共数据模型和输出逻辑。

#### 验收标准

- 拆分前后的 ID、BED、position TSV 和 SNP calling BED 字节级一致。
- 现有单元测试和端到端测试全部通过。
- 每个模块有单一职责，CLI 层不再包含核心 DataFrame 变换逻辑。

## 推荐实施顺序

1. 修复 Conda PATH，并用当前番茄 fixture 恢复 ID/VCF smoke test。
2. 固化依赖环境，将 smoke test 接入 CI。
3. 完成帮助、schema 和参数边界校验。
4. 增加已有 chain 复用入口和 chain 元数据。
5. 建立大型基因组基线，再实施 PAF、资源和窗口 overlap 优化。
6. 在行为测试稳定后拆分 Python 辅助脚本。

P0 是后续所有性能测试和功能验收的前置条件。P1 项目宜在继续扩展流程功能前完成；P2/P3 应以实际基准和维护成本为依据，不建议一次性重写。

## 建议拆分的 Issue

| 建议标题 | 优先级 | 依赖 |
| --- | --- | --- |
| 修复 Nextflow process 环境被 Miniforge base PATH 覆盖 | P0 | 无 |
| 增加锁定版本的 Conda/Apptainer 运行环境 | P1 | P0 |
| 为 ID 与 VCF 模式建立 CI smoke test | P1 | P0、可复现环境 |
| 实现 `--help` 与统一参数 schema 校验 | P1 | 无 |
| 支持通过 `--chain` 跳过基因组比对 | P1 | P0、端到端测试 |
| 增加可选 PAF 发布与资源基准 | P2 | CI、代表性数据 |
| 评估并实现 split window overlap | P2 | 边界基准数据 |
| 修复重复运行时 report/timeline 覆盖策略 | P2 | 无 |
| 拆分 liftover Python 数据处理与 CLI | P3 | 端到端测试 |

## 验证记录

### 评估时（2026-07-15）

```text
Python unittest: 21 tests, OK
Nextflow config -profile standard: exit 0
Missing required input: exit 1, correctly reports --ref_fa
ID smoke: reaches LIFTOVER_BY_ID, then fails because base Python shadows task environment
```

### 实施后（2026-07-16）

按推荐顺序已落地：P0 PATH 修复；P1 可复现环境、CI、help/参数校验、`--chain` 复用；P2 report 覆盖、`publish_paf`、large_mem 重试增长。P2 滑窗 overlap 与 P3 Python 拆分留待基准与行为锁定后进行。

审核反馈修复（同日）：

- `WRITE_CHAIN_META` 改为 base64/JSON + 引用 heredoc；校验逻辑抽到 `bin/chain_meta.py`
- reuse 校验 chain header 染色体名/长度/方向，可选 `--chain_meta` sidecar sha256
- reuse 发布 `chain/all.chain` 副本
- CI：`cp nextflow.config.example nextflow.config` + `-profile ci`；严格 `*_path` 断言；trace 证明未跑 alignment
- `environment.yml` 固定版本；`run_meta.yml` 记录 git commit/dirty
- 全模块 versions 片段改用 echo 写入；collate 过滤字面量 `END`

```text
Python unittest: 26 tests, OK
nextflow run . --help (clean checkout + example config): exit 0
ID smoke: success; software_versions.yml valid (no END); tool paths under nf-liftover-tools
ID + --chain + --chain_meta: success; publishes chain/all.chain; no ALIGN_* in trace
Swapped ref/query with same chain: fails with "Chain appears reversed"
run_meta: git_commit + git_dirty recorded
```
