# `--aligner` 开关收尾 Spec

## 背景

`docs/wheat-split-benchmark-spec.md` 的阶段 D 已经给出结论：minimap2 与 mm2plus 不存在普适最优，交叉点在并发投递数 P≈8-12，因此以 `--aligner` 开关交给用户按执行环境选择，默认保持 `minimap2`。开关已在工作树中实现（8 个文件，+53 -12 行）：新增 `tool_aligner` label、`params.aligner` / `params.aligner_envs`、三个比对 process 的 aligner 内插，以及 `main.nf` 的取值校验。

功能本身可用，但收尾有缺口。最主要的一条是：`nextflow.config` 是 gitignored 的用户私有文件，新参数只落在 `nextflow.config.example`，存量用户拉到该改动后会在启动阶段直接崩溃。本 spec 约束合并前必须完成的修复与收尾项。

## 目标

1. 存量 `nextflow.config`（无 `aligner` / `aligner_envs`）拉到本改动后仍能按原行为运行，不需要用户手工改配置。
2. 用户显式请求了配置无法提供的 aligner 时，报错可操作，能指向修复方式。
3. README 记录 `--aligner`，使新增的用户可见开关有文档。
4. `--aligner` 的两个取值与非法取值都有自动化测试覆盖。
5. 纠正 `docs/mm2plus-benchmark-spec.md` 中已被后续实验推翻的非目标声明。

## 非目标

- 不改变默认行为：默认仍是 `minimap2`，默认路径的输出必须逐字节不变。
- 不调整 `split_mem` 资源（`cpus` / `memory`）。基准文档指出这是收益更大的独立议题（P=16 的 2.59 h vs 当前 P=4 的 6.41 h，约 2.5x），但应单开一轮评估，不混入本次。
- 不接入 `nf-schema` 插件（见 P2 观察项，属独立议题）。
- 不新增第三个 aligner，不改动 `conf/modules.config` 的 `ext.args`。
- 不重跑基准。阶段 A-D 数据已定稿，本轮只做工程收尾。

## 问题清单

### P0: 存量 `nextflow.config` 缺少新参数导致启动即崩

`.gitignore:14` 将 `/nextflow.config` 排除在版本控制外，用户配置由 `nextflow.config.example` 拷贝而来且此后独立演进。本次新增的两个参数只写在 example（`nextflow.config.example:18` 和 `:34-37`），存量用户的配置里不存在。

证据：

- `main.nf:23` 执行 `params.aligner_envs.containsKey(params.aligner)`。`params.aligner_envs` 为 null 时抛 `Cannot invoke method containsKey() on null object`，报错既不指向参数名也不指向修复方式。
- `conf/conda.config:12` 的 `beforeScript` 内插 `params.aligner_envs[params.aligner]`，同样在 null 上取下标。
- 三个比对 module（`align_split_window.nf:30`、`align_whole_chromosome.nf:39`、`build_query_mmi.nf:35`）的 `def aligner = params.aligner` 会得到 null，脚本退化为 `null -cx asm5 ...`。
- `nextflow_schema.json:49-53` 虽然给了 `"default": "minimap2"`，但该 schema 当前不生效（见 P2），无法兜底。

要求（方案 A，优先）：把兜底默认值下沉到仓库跟踪的 `conf/base.config`，用自引用写法避免覆盖用户已有设置。`nextflow.config.example:40` 的 `includeConfig 'conf/base.config'` 位于 params 块之后，因此 base.config 中的直接赋值会**覆盖**用户配置，必须写成回退形式：

```groovy
// conf/base.config
params {
    // 存量 nextflow.config 可能没有这两项；此处仅在未定义时回退，
    // 不覆盖用户在自己 params 块中的显式设置。
    aligner      = params.aligner      ?: 'minimap2'
    aligner_envs = params.aligner_envs ?: [ minimap2: params.conda_envs?.ngs ?: 'nf-liftover-tools' ]
}
```

该写法必须先验证再采纳，验证方式见下方验收。若实测顺序语义与预期不符（用户显式设置被覆盖，或回退未生效），改用方案 B：在每个取值点就地回退，不动 base.config。

- `main.nf:23`：改为 `def alignerEnvs = params.aligner_envs ?: [:]`、`def aligner = params.aligner ?: 'minimap2'`，仅当 `aligner != 'minimap2' && !alignerEnvs.containsKey(aligner)` 时抛错。
- `conf/conda.config:12`：改为 `params.aligner_envs?.get(params.aligner ?: 'minimap2') ?: params.conda_envs.ngs`。
- 三个 module：`def aligner = params.aligner ?: 'minimap2'`（script 与 stub 各一处，共 6 处）。

无论采用哪个方案，`main.nf` 的错误文案都要可操作，至少包含受支持取值与配置来源：

```
Unsupported --aligner 'xxx'; supported: minimap2, mm2plus.
If you upgraded an existing nextflow.config, copy the `aligner_envs` block from nextflow.config.example.
```

验收：

- 构造一份不含 `aligner` / `aligner_envs` 的配置（可由 `nextflow.config.example` 删去这两项得到），以 `-stub-run` 跑通 tomato smoke，退出码 0，且 `ALIGN_*` 任务的 `.command.sh` 中命令为 `minimap2`。
- 同一份旧配置加 `--aligner mm2plus`，报错信息包含受支持取值列表与 `nextflow.config.example` 字样，不出现 `NullPointerException`。
- 采用方案 A 时另需验证不覆盖：在 `nextflow.config` 的 params 块里写 `aligner = 'mm2plus'`，不加任何 CLI 参数，`-stub-run` 后 `.command.sh` 中为 `mm2plus`；CLI `--aligner minimap2` 能进一步覆盖为 `minimap2`。

### P1: README 未记录 `--aligner`

新增的是用户可见开关，README 是唯一的用户入口文档，当前零提及。

证据：

- `README.md:17` 的流程图把 `minimap2` 写死在链路上。
- `README.md:28` 的运行要求列出 `minimap2`，未说明可选 mm2plus 及其额外环境。

要求：

- 流程图中 `minimap2` 改为中性表述（如 `比对`），或标注为 `minimap2/mm2plus`。
- 运行要求中说明 mm2plus 为可选，需要额外的 conda 环境，环境名由 `params.aligner_envs` 指定。
- `README.md:103` 的 `## 主要参数` 表格中，在 `--pair_strategy` 行之后插入 `--aligner` 行（与 `nextflow.config.example` 的 params 顺序一致），默认值 `minimap2`，说明列写明可选 `mm2plus`。
- 表格下方补一句选型依据：窗口数远多于可用并发槽位时选 `minimap2`，窗口数少于可用槽位（少染色体 liftover、独占大节点）时选 `mm2plus`，细节链接到 `docs/wheat-split-benchmark-spec.md` 的最终结论表。

验收：

- `grep -n 'aligner' README.md` 有命中，且包含默认值与两个可选值。
- README 中不再存在把 minimap2 表述为唯一比对器的句子。

### P1: `--aligner` 无测试覆盖

`nf-test.config` 的 `configFile` 指向 `nextflow.config.example`、`profile` 为 `test`，因此新参数在测试环境中天然可见；`conf/test.config:30-34` 也已经按 `params.aligner` 准备好两套 bioconda 依赖。当前 `tests/tomato_smoke.nf.test` 的所有 case 都跑默认路径，mm2plus 分支与非法取值分支均未被执行过。

要求：

- 新增一个 mm2plus smoke case，参数与既有「runs tomato fixture end to end」一致，仅追加 `aligner = "mm2plus"`。放进独立文件 `tests/tomato_smoke_mm2plus.nf.test`，使本机未安装 mm2plus 时可以只跑 `nf-test test tests/tomato_smoke.nf.test` 而不受影响（仓库当前无 CI 配置，测试由人工本地执行，故不引入 tag 或 profile 机制）。
- 在 `tests/tomato_smoke.nf.test` 中新增一个非法取值 case，如 `aligner = "bwa"`，断言 `workflow.failed` 且 stdout 含受支持取值列表。该 case 不依赖 mm2plus，留在主文件。
- mm2plus case 的断言**不得**对 `chain/all.chain` 做逐字节或哈希比较。`docs/mm2plus-benchmark-spec.md:235` 已指出 chain 文件跨 aligner 不可 byte 比较（chain id 编号与记录顺序会变）。断言限于：`workflow.success`、`trace.failed().size() == 0`、`all.chain` 存在且非空，以及 liftover 产物（`.id` / `.bed` / `.snpcalling.bed`）文本与默认 case 相同——后者是真正的等价性判据，坐标级结果应当一致。
验收：

- `nf-test test tests/tomato_smoke.nf.test` 全绿，用例数比当前增加 1（非法取值 case）。
- `nf-test test tests/tomato_smoke_mm2plus.nf.test` 全绿（需本机可安装 mm2plus）。
- mm2plus case 与默认 case 的三个 liftover 产物文本断言完全一致。
- 移除 `main.nf` 的 aligner 校验分支后，非法取值 case 变红（反向验证该用例确实在测这段逻辑）。

### P1: `mm2plus-benchmark-spec.md` 的非目标已被推翻

`docs/mm2plus-benchmark-spec.md:24` 写着「不新增 `--aligner` 开关」。该判断基于水稻单条染色体数据（1.20x，结论倾向不接入），后续小麦切窗口实验推翻了它并落地了开关。两份文档现在互相矛盾，按时间顺序读会误导。

要求：

- 在 `docs/mm2plus-benchmark-spec.md` 的非目标条目或结论区加一行 superseded 注记，说明该条仅约束当轮实验范围，最终决策见 `docs/wheat-split-benchmark-spec.md` 的最终结论。
- 不修改该文档已记录的实验数据与当轮结论，只加注记。

验收：

- `docs/mm2plus-benchmark-spec.md` 中出现指向 `wheat-split-benchmark-spec.md` 的交叉引用。
- 两份文档不再出现无限定语的「不新增 `--aligner` 开关」表述。

### P2（观察项，本轮不实施）: `nextflow_schema.json` 未接入

全仓库无 `plugins` 块，`main.nf` 也没有 `validateParameters()` 调用：

```
grep -rn "nf-schema\|validateParameters\|plugins" --include='*.nf' --include='*.config' --include='*.example' .   # 无输出
```

因此 `nextflow_schema.json` 目前是纯文档，其中的 `"default"` 与 `"enum"` 都不生效——本次给 `aligner` 写的 `"default": "minimap2"` 和 `"enum": ["minimap2","mm2plus"]` 一条都没有运行时作用，参数校验只能靠 `main.nf` 手写。这也是 P0 需要额外兜底的直接原因。

`docs/remediation-spec.md` 曾以「项目已经使用 nf-schema」为前提删除手写 help，与当前事实不符，该前提应在接入 nf-schema 时一并复核。

本轮不实施，仅记录。接入后可移除 `main.nf:23-26` 的手写校验，由 schema enum 接管。

## 已核对、无需改动

以下几处在审核中确认不受 `tool_ngs` → `tool_aligner` label 变更影响，列出以免后续重复排查：

- `conf/slurm.config` 只设置 executor / queue 与 `executor` 块，无任何 `withLabel` 覆盖，label 改名不丢配置。
- `conf/modules.config` 的 `ext.args`（`-x asm5`、`-cx asm5 --cs`）对两个 aligner 通用，mm2plus 为 drop-in fork。
- `conf/test.config:31` 的 `params.aligner == 'mm2plus' ? ... : ...` 对 null 安全，未定义时落到 minimap2 分支，与 P0 的回退方向一致。
- `-resume` 跨 aligner 切换不会误用缓存：`params.aligner` 被内插进 process script，取值变化即改变任务哈希，`BUILD_QUERY_MMI` 与 `ALIGN_*` 会正确重算，不存在 minimap2 建的 `.mmi` 喂给 mm2plus 的情况。
- 三个 module 的 stub 由 `minimap2: 2.26` 改为 `${aligner}: stub`，现有测试只断言 `software_versions.yml` 存在，不受影响。

## 实施阶段

### 阶段 1: P0 兼容性修复

1. 先验证方案 A 的配置顺序语义（用一份显式设 `aligner = 'mm2plus'` 的配置 `-stub-run` 观察 `.command.sh`）。
2. 按验证结果落地方案 A 或方案 B。
3. 改写 `main.nf` 错误文案。
4. 跑通 P0 的三条验收。

### 阶段 2: 测试覆盖

5. 在 `tests/tomato_smoke.nf.test` 补非法取值 case，新建 `tests/tomato_smoke_mm2plus.nf.test` 补 mm2plus case。
6. 执行反向验证（临时移除 `main.nf` 的校验分支，确认非法取值用例变红，随后还原）。

### 阶段 3: 文档

7. 更新 README。
8. 给 `mm2plus-benchmark-spec.md` 加 superseded 注记。

### 阶段 4: 提交

9. 拆两个 commit：实现（`main.nf`、三个 module、四个 config、`nextflow_schema.json`、README、测试）与基准文档（`docs/mm2plus-benchmark-spec.md` 注记、`docs/wheat-split-benchmark-spec.md`、本 spec）。后者体量大且是独立可引用的证据，混在一起会淹没实现改动的 review。

## 最终验收标准

- 存量配置（无新参数）`-stub-run` 通过，命令为 `minimap2`，无 `NullPointerException`。
- 存量配置 + `--aligner mm2plus` 报错可操作，含受支持取值与配置来源。
- 两个测试文件均全绿，共新增 2 个用例；mm2plus 与 minimap2 的 liftover 产物文本一致。
- 默认路径产物与本次改动前逐字节一致（`chain/all.chain` 与三个 liftover 文件）。
- README 记录 `--aligner`，含默认值、可选值与选型依据链接。
- 两份基准文档不再互相矛盾。

## 风险与备注

- 方案 A 依赖 Nextflow 对跨文件 params 自引用的求值顺序。这是本 spec 中唯一未经实测的假设，必须先验证再落地，否则会静默覆盖用户配置——这种失败不报错，只是行为不对，比崩溃更难发现。
- mm2plus 测试用例引入了对 bioconda `mm2plus>=1.3` 的网络依赖，首次运行需要现建 conda 环境，会显著变慢。这是把它拆成独立文件的原因，不要为了省时间删掉该用例。
- 阶段 D 的交叉点（P≈8-12）来自单机 32 核 / 122 GB 环境。README 的选型建议应表述为「取决于并发投递数与可用槽位的关系」，而不是复制具体 P 值，避免用户在不同规格的机器上照搬。
- `split_mem` 右调按用户要求维持现状。若后续调整 `cpus`，会直接移动并发投递数，进而移动 aligner 的最优选择——两个议题的耦合关系需要在那一轮的 spec 里显式处理。
