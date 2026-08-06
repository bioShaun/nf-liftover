# Drift Audit Report — 2026-08-05

文档与代码/配置/现实之间的偏差清单。每一条都经过 Read/Grep/`test -f` 对照仓库当前状态核实，行号已重新校验。

## Summary

- **31** 个 drift item，分布在 **10** 个文件中。
- Action 分布：mark as legacy-known-drift (no fix) × 8；annotate as superseded/planned × 11；fix doc × 7；no action (historical benchmark artifact) × 2；external repo, verify separately × 1；fix doc (env name) × 1；fix doc (line refs) × 1。
- Scope：全部为 in-repo（本次只修仓库内文档，外部仓库只读）。
- 审计基准：worktree `worktree-rapid-harbor-dd26`，2026-08-05。

## Follow-up 2026-08-06：NF26 迁移落地后的状态更新

审计完成后，Nextflow 24→26 迁移收尾于同日落地，以下审计事实与条目状态随之改变（保留上文原文供历史追溯）：

- 「全仓库无 `plugins` 块」不再成立：`nextflow.config.example` 与 `nextflow.legacy.config` 均声明 `plugins { id 'nf-schema@2.7.2' }`，`main.nf` 已 `include ... from 'plugin/nf-schema'`，并以 `schemaValidate()` 作手写校验之后的 schema 级第二道门。**无 try/catch 静默降级**：插件由 config `plugins` 块声明，加载失败会在启动解析阶段中止运行，`main.nf` 内 catch 不到也不应 catch（2026-08-06 复核时发现原实现的 try/catch 会同时吞掉真实校验错误，已移除）。
- 受影响条目：#11/#18/#22 的「未接入」判定仅对审计时点成立；`docs/remediation-spec.md` 勘误与 `docs/usage.md:149` 的表述已于 2026-08-06 相应更新。#11 grep 命中的 `upgrade-plan.md:358`「已接入 nf-schema」自此成为真实陈述，条目关闭。
- publishDir 已全部迁移至 workflow `output:` 块（7 模块删 8 处），`upgrade-plan.md` L124/L197-207/L329/L359 的 output 块规划按 NF26 target-based 语法落地；#17 条目随之关闭。
- 下方 Post-Fix Re-Grep Checks 中 #11/#18/#22/#23/#24/#31 的期望已按实际修复策略（annotate / 保留原文）于 2026-08-06 修正——原「期望：无输出」与保留原文策略矛盾，且 #24 正则会误匹配已更正的 `tc-probe-design-dev`。
- 验收（2026-08-06 复核重跑）：39 python 单测、12/12 nf-test、3 profile config parse、-stub 双向发布（`publish_paf=false/true`）全部通过。

## Verified Machine Facts

以下事实在审计时通过 `ls`/`test`/`mamba run` 实地核对，作为 drift 判定的 ground truth。

### Nextflow profiles

- `nextflow.config.example:64-88` 定义 5 个 profile：`standard`、`slurm`、`slurm_new`、`test`、`ci`。
- `nextflow.legacy.config:53-71` 定义 4 个，除 `ci` 外与其余一致。
- `slurm` → queue `normal`，`slurm_new` → queue `cae`（`conf/slurm.config` + `nextflow.config.example:73-76`）。
- `slurm_new2`：在所有 `.config` / `.example` / `.nf` 文件中零匹配，仅出现在 `docs/pipe.md:85,97` 和 `docs/upgrade-plan.md:18,174`。

### Conda 环境

- `/project/software/miniforge3/envs/` 下存在：`nextflow`、`nextflow26`、`nf-liftover-tools`、`bioinfo`、`py-13`、`mm2plus` 等。
- **不存在**：`nextflow24`、`ngs`、`py3-13`、`probe-design`。
- `tc-probe-design-dev` 存在（`probe-design` 不存在）。
- `nf-liftover-tools` 环境内容（`environment.yml`）：`python`、`pandas`、`typer`、`loguru`、`pyfaidx`、`samtools`、`htslib`、`minimap2`、`transanno`。**无 `seqkit`**（`mamba run -n nf-liftover-tools which seqkit` 无输出）。
- `nextflow26` 环境中 `nf-test 0.9.5` 已安装（`mamba run -n nextflow26 nf-test version` 输出 `nf-test 0.9.5`）。

### 仓库结构

- `modules/local/` 下 16 个 `.nf` 文件，**无** `restore_split_paf.nf`，有 `combine_split_pafs.nf`。
- `assets/` 目录不存在。
- `bin/liftover_by_id.py` 和 `bin/merge_blast_rescue.py` 均存在。
- `examples/sl4-vs-la2093/` 仅含 `run.sh`、`README.md`、`chrom_pairs.tsv`，**无** `TCZZSL20K.id`。
- `.github/workflows/ci.yml` 存在。
- `main.nf` 有手写 `validateParameters()`（:163）和 `helpMessage()`（:13），help 分支用 `return`（:262-265），**无** `System.exit(0)`。全仓库无 `plugins` 块。
- `conf/base.config:45`：`large_mem` memory = `{ 20.GB * task.attempt }`（非 60.GB）。
- `nextflow.config.example:36-39`：`conda_envs` 映射 `ngs→nf-liftover-tools`、`py→nf-liftover-tools`（非 `bioinfo`/`py-13`）。
- `nextflow.config.example:41-44`：`aligner_envs` 映射 `minimap2→nf-liftover-tools`、`mm2plus→mm2plus`。
- `nf-test.config:6`：`configFile "nextflow.legacy.config"`（非 `nextflow.config.example`）。
- `docs/usage.md` 全文无 `nf-schema` 插件安装步骤（:149 明确写"不依赖 `nf-schema` 插件"）。

### 文件系统路径

| 路径 | 状态 |
|------|------|
| `/public/scripts/nf-liftover` | 存在 |
| `/public/scripts/nf-liftover/examples/sl4-vs-la2093/chrom_pairs.tsv` | 存在 |
| `/public/data/genomes`（复数） | 不存在 |
| `/public/data/genomes/solanum_lycopersicum_LA2093` | 不存在 |
| `/public/data/genome/solanum_lycopersicum_LA2093`（单数） | 不存在 |
| `/project/tmp/mm2plus-bench/` | 不存在 |
| `/project/tmp/wheat-split-bench/` | 不存在 |
| `/public/home/zxchen/software/miniconda3` | 不存在 |
| `~/scripts/ngs-utils-v24` | 不存在 |
| `~/scripts/tc-pytools` | 不存在 |
| `~/scripts/generic/probe-design` | 不存在 |
| `/data_0/liftover/arachis_hypogaea_Fuhuasheng-to-arachis_hypogaea_tifrunner` | 存在 |
| `/data_0/panel_design/projects/TC-BY-peanut` | 存在 |
| `~/.nextflow/plugins/nf-schema-2.7.2` | 存在（upgrade-plan 目标 2.2.0） |

### Issue 文档验证

- `docs/issue-001/002/003` 中引用的 `/data_0/` 路径均已验证存在。
- `bin/liftover_by_id.py` + `bin/merge_blast_rescue.py` 存在。
- `tc-pytools/panel/realign_blast.py` 属外部仓库，本机不存在，无法验证（external, verify separately）。

## Drift Items

| # | File:Line | Claim | Reality | Action | Scope |
|---|-----------|-------|---------|--------|-------|
| 1 | docs/pipe.md:85 | `-profile slurm_new2` | `slurm_new2` 在所有 config/code 中零匹配；只有 `slurm`/`slurm_new` | mark as legacy-known-drift (no fix) | in-repo |
| 2 | docs/pipe.md:97 | `-profile slurm_new2` | 同 #1 | mark as legacy-known-drift (no fix) | in-repo |
| 3 | docs/pipe.md:28 | `conda activate ngs` | conda env `ngs` 不存在 | mark as legacy-known-drift (no fix) | in-repo |
| 4 | docs/pipe.md:74 | `conda activate nextflow24` | conda env `nextflow24` 不存在 | mark as legacy-known-drift (no fix) | in-repo |
| 5 | docs/pipe.md:73 | `source /public/home/zxchen/software/miniconda3/...` | 路径不存在 | mark as legacy-known-drift (no fix) | in-repo |
| 6 | docs/pipe.md:115,151 | `conda activate py3-13` | 实际 env 名为 `py-13`（无 `3`） | mark as legacy-known-drift (no fix) | in-repo |
| 7 | docs/pipe.md:123,134 | `~/scripts/tc-pytools/liftover/liftover_by_id.py` | 路径不存在；vendored 副本在 `bin/liftover_by_id.py` | mark as legacy-known-drift (no fix) | in-repo |
| 8 | docs/pipe.md:80-86,92-98 | `~/scripts/ngs-utils-v24/workflows/align_chromosomes*.nf` + `--chr_pairs`/`--species_a_dir`/`--species_b_dir` flags | 路径不存在；flags 属旧工作流，当前流水线无此参数 | mark as legacy-known-drift (no fix) | in-repo |
| 9 | docs/pipe.md:153-163 | `~/scripts/generic/probe-design/prepare-probe-dir.sh` | 路径不存在；README:252 已声明不属于本项目 | mark as legacy-known-drift (no fix) | in-repo |
| 10 | docs/upgrade-plan.md:18 | `slurm_new2` profile 偏差 | 自述已知偏差；`slurm_new2` 确实不存在 | annotate as superseded/planned | in-repo |
| 11 | docs/upgrade-plan.md:87,120,356 | nf-schema `validateParameters`/`paramsSummaryLog` "已接入" | :356 标记 `[x]` 已完成，但全仓库无 `plugins` 块；`main.nf:163` 是手写 `validateParameters()`。:87/:120 在目标规划段描述 target state，非完成声明 | annotate as superseded/planned | in-repo |
| 12 | docs/upgrade-plan.md:91,224-225 | `conda_envs` 映射到 `bioinfo`/`py-13` | 实际 `nextflow.config.example:36-39`：`ngs→nf-liftover-tools`、`py→nf-liftover-tools` | annotate as superseded/planned | in-repo |
| 13 | docs/upgrade-plan.md:98-101 | 目标 modules `seqkit/split.nf`、`seqkit/sliding.nf` 等 | 未落地；`modules/local/` 下 16 个文件均为 local/ 命名，无 seqkit/ 前缀 | annotate as superseded/planned | in-repo |
| 14 | docs/upgrade-plan.md:102 | `local/restore_split_paf.nf` | 不存在；等价逻辑在 `modules/local/combine_split_pafs.nf` | annotate as superseded/planned | in-repo |
| 15 | docs/upgrade-plan.md:107-108 | `assets/` 目录含 `nextflow_schema.json`、`schema_input.json` | `assets/` 目录不存在；`nextflow_schema.json` 在仓库根 | annotate as superseded/planned | in-repo |
| 16 | docs/upgrade-plan.md:161 | `--by-name-suffix` flag | 实际 CLI 为 `--pair_strategy {order,suffix}`（`nextflow.config.example:20`，`main.nf:200-203`） | annotate as superseded/planned | in-repo |
| 17 | docs/upgrade-plan.md:197-207 | workflow `output:` 块语法 `directory "..." { 'x' { from ... } }` | NF26 24-04 first-preview 形式已过时；当前用 `publishDir`（:357 自述未切换） | annotate as superseded/planned | in-repo |
| 18 | docs/upgrade-plan.md:248,363 | `docs/usage.md` 写明 `nextflow plugin install nf-schema@2.2.0` | `docs/usage.md` 全文无此内容；:149 明确写"不依赖 `nf-schema` 插件" | annotate as superseded/planned | in-repo |
| 19 | docs/upgrade-plan.md:285-290,335 | `/public/data/genomes/solanum_lycopersicum_LA2093/` "已就位" | 目录不存在（复数和单数均无） | annotate as superseded/planned | in-repo |
| 20 | docs/upgrade-plan.md:313 | `examples/sl4-vs-la2093/TCZZSL20K.id` 软链 | 文件不存在；目录仅含 `run.sh`、`README.md`、`chrom_pairs.tsv` | annotate as superseded/planned | in-repo |
| 21 | README.md:29 | 任务环境含 `seqkit` | `environment.yml` 和 `nf-liftover-tools` env 均无 seqkit | fix doc | in-repo |
| 22 | docs/remediation-spec.md:5,12,51,55 | "项目已经使用 `nf-schema`" + 删除 `System.exit(0)` help 分支 | nf-schema 未接入；`main.nf` help 分支用 `return`（:262-265），无 `System.exit(0)` | annotate as superseded/planned | in-repo |
| 23 | docs/remediation-spec.md:86 | `large_mem` 使用 `{ 60.GB * task.attempt }`，第三次 180GB | 实际 `conf/base.config:45`：`{ 20.GB * task.attempt }`（20→40→60） | annotate as superseded/planned | in-repo |
| 24 | docs/remediation-spec.md:281,297 | `mamba run -n probe-design` | env `probe-design` 不存在；机器上有 `tc-probe-design-dev` | annotate as superseded/planned | in-repo |
| 25 | docs/aligner-switch-followup-spec.md:92 | `nf-test.config` 的 `configFile` 指向 `nextflow.config.example` | 实际 `nf-test.config:6`：`configFile "nextflow.legacy.config"` | fix doc | in-repo |
| 26 | docs/aligner-switch-followup-spec.md:96 | "仓库当前无 CI 配置" | `.github/workflows/ci.yml` 存在 | fix doc | in-repo |
| 27 | docs/aligner-switch-followup-spec.md:29 | 行号引用 `nextflow.config.example:18` 和 `:34-37` | 实际 `aligner` 在 :21，`aligner_envs` 在 :41-44（偏移 3-7 行） | fix doc | in-repo |
| 28 | docs/optimization-assessment.md:17 | "`nextflow26` 环境中没有可直接调用的 `nf-test`" | 已过时；`nf-test 0.9.5` 现已安装在 `nextflow26` env | fix doc | in-repo |
| 29 | docs/mm2plus-benchmark-spec.md:26,46-52,71,77,108,118,180 | `/project/tmp/mm2plus-bench/` 产物路径 | 目录已不存在（:41 正确指出番茄数据目录不存在） | no action (historical benchmark artifact) | in-repo |
| 30 | docs/wheat-split-benchmark-spec.md:51,63-79,193,268-272 | `/project/tmp/wheat-split-bench/` 产物路径 | 目录已不存在 | no action (historical benchmark artifact) | in-repo |
| 31 | examples/sl4-vs-la2093/run.sh:10-12, examples/sl4-vs-la2093/README.md:5-7 | `/public/data/genomes/solanum_lycopersicum_LA2093/` 下三个文件 | 目录不存在（复数和单数均无）；`run.sh:13` 的 mapping 路径在部署副本中存在 | fix doc | in-repo |

### 附加发现（env 名称不一致）

| # | File:Line | Claim | Reality | Action | Scope |
|---|-----------|-------|---------|--------|-------|
| 32 | examples/sl4-vs-la2093/run.sh:6, examples/sl4-vs-la2093/README.md:12 | `conda activate nextflow` / `envs/nextflow/bin/nextflow` | README.md:57、docs/usage.md:13、tests/data/tomato-smoke/run-smoke.sh:12 均用 `nextflow26`；两个 env 都存在但名称不一致 | fix doc | in-repo |

### 非 drift 的已验证项

以下在 inventory 中提及但经核实无 drift，记录在此供后续审计跳过：

- `docs/issue-001/002/003` 中 `/data_0/` 路径：已验证存在。
- `bin/liftover_by_id.py`、`bin/merge_blast_rescue.py`：存在。
- `tc-pytools/panel/realign_blast.py`：外部仓库，本机不存在，无法验证。Action: external repo, verify separately。Scope: in-repo（本仓库无 drift）。
- `examples/sl4-vs-la2093/run.sh:13` mapping 路径 `/public/scripts/nf-liftover/examples/sl4-vs-la2093/chrom_pairs.tsv`：在部署副本中存在。

## Post-Fix Re-Grep Checks

修复完成后运行以下命令验证 drift 已消除。每条标注对应的 item 编号。

```bash
# --- #1,#2: slurm_new2 不应出现在 pipe.md（保留 upgrade-plan 自述则跳过） ---
grep -n 'slurm_new2' docs/pipe.md
# 期望：无输出（或标注 legacy 后保留但加注记）

# --- #3: ngs env 不应出现在 pipe.md ---
grep -n 'conda activate ngs' docs/pipe.md
# 期望：无输出

# --- #4: nextflow24 不应出现在 pipe.md ---
grep -n 'nextflow24' docs/pipe.md
# 期望：无输出

# --- #5: /public/home/zxchen 不应出现在 pipe.md ---
grep -n '/public/home/zxchen' docs/pipe.md
# 期望：无输出

# --- #6: py3-13 不应出现在 pipe.md ---
grep -n 'py3-13' docs/pipe.md
# 期望：无输出

# --- #7: ~/scripts/tc-pytools 不应出现在 pipe.md ---
grep -n '~/scripts/tc-pytools' docs/pipe.md
# 期望：无输出

# --- #8: ngs-utils-v24 不应出现在 pipe.md ---
grep -n 'ngs-utils-v24' docs/pipe.md
# 期望：无输出

# --- #9: prepare-probe-dir.sh 在 pipe.md 中应标注 legacy/out-of-scope ---
grep -n 'prepare-probe-dir' docs/pipe.md
# 期望：保留但带 legacy 注记

# --- #11: nf-schema "已接入" 应改为 "计划中" 或删除 [x] ---
grep -n 'nf-schema.*已接入\|schema-polish.*已接入' docs/upgrade-plan.md
# 期望：无输出
# 2026-08-06 修正：nf-schema 已真实接入，命中 :358 现为真实陈述，条目关闭（见文首 Follow-up）

# --- #12: conda_envs 不应映射到 bioinfo/py-13 in upgrade-plan ---
grep -n "ngs : 'bioinfo'\|py  : 'py-13'" docs/upgrade-plan.md
# 期望：无输出（或标注为历史规划）

# --- #14: restore_split_paf.nf 不应作为目标模块出现 ---
grep -n 'restore_split_paf\.nf' docs/upgrade-plan.md
# 期望：无输出或标注 superseded

# --- #15: assets/ 目录不应在 upgrade-plan 中作为目标结构 ---
grep -n '^- `assets/`' docs/upgrade-plan.md
# 期望：无输出或标注 superseded

# --- #16: --by-name-suffix 不应出现 ---
grep -n 'by-name-suffix' docs/upgrade-plan.md
# 期望：无输出

# --- #18: usage.md 不应被声称含 nf-schema 安装步骤 ---
grep -n 'usage.md.*plugin install\|usage.md.*nf-schema' docs/upgrade-plan.md
# 期望：无输出（命中处为历史规划段，属"保留原文"策略可接受）
# 2026-08-06 修正：usage.md 现已如实描述 nf-schema 为启动依赖及离线缓存要求，原 drift 方向已反转

# --- #19: /public/data/genomes "已就位" 应修正 ---
grep -n '已就位\|已就位.*genomes' docs/upgrade-plan.md
# 期望：无输出或标注路径不存在

# --- #20: TCZZSL20K.id 软链声明应修正 ---
grep -n 'TCZZSL20K.id.*软链' docs/upgrade-plan.md
# 期望：无输出或标注文件不存在

# --- #21: README 不应声称 seqkit 在任务环境中 ---
grep -n 'seqkit' README.md
# 期望：无输出

# --- #22: remediation-spec 不应声称 nf-schema 已使用 ---
grep -n '已经使用.*nf-schema\|项目已经使用' docs/remediation-spec.md
# 期望（2026-08-06 修正）：命中仅限保留的原文（:54），且同节均带 2026-08-05/06 勘误标注

# --- #22: remediation-spec 不应声称 System.exit(0) 存在 ---
grep -n 'System.exit(0)' docs/remediation-spec.md
# 期望（2026-08-06 修正）：命中仅限保留的原文/目标引用（:14/:54/:295），且均带勘误标注

# --- #23: large_mem 不应写 60.GB ---
grep -n '60\.GB.*task.attempt' docs/remediation-spec.md
# 期望（2026-08-06 修正）：命中仅限保留的原文（:91），其后紧随 2026-08-05 勘误段

# --- #24: probe-design env 不应出现在 remediation-spec ---
grep -n 'probe-design' docs/remediation-spec.md | grep -v 'tc-probe-design-dev'
# 期望（2026-08-06 修正）：无输出（原正则会误匹配已更正的 tc-probe-design-dev，需排除）

# --- #25: nf-test.config 不应被声称指向 nextflow.config.example ---
grep -n 'nextflow.config.example' docs/aligner-switch-followup-spec.md
# 期望：无输出（或改为 nextflow.legacy.config）

# --- #26: "仓库当前无 CI 配置" 不应出现 ---
grep -n '无 CI 配置' docs/aligner-switch-followup-spec.md
# 期望：无输出

# --- #27: 行号引用应更新 ---
grep -n 'nextflow.config.example:18\|nextflow.config.example:34-37' docs/aligner-switch-followup-spec.md
# 期望：无输出（改为 :21 和 :41-44）

# --- #28: "nextflow26 环境中没有 nf-test" 不应出现 ---
grep -n '没有.*nf-test\|中没有.*nf-test' docs/optimization-assessment.md
# 期望：无输出或标注已过时

# --- #29,#30: /project/tmp benchmark 目录（历史产物，不修） ---
test -d /project/tmp/mm2plus-bench && echo "EXISTS" || echo "GONE"
test -d /project/tmp/wheat-split-bench && echo "EXISTS" || echo "GONE"
# 期望：GONE / GONE（no action，历史基准产物）

# --- #31: /public/data/genomes 路径不应在 examples 中出现 ---
grep -rn '/public/data/genomes' examples/sl4-vs-la2093/
# 期望（2026-08-06 修正）：命中仅限保留的示例路径，且 README:9 与 run.sh:9-19 均带 2026-08-05 审计警告注记

# --- #32: env 名称应统一为 nextflow26 ---
grep -rn 'conda activate nextflow\b' examples/sl4-vs-la2093/
grep -rn 'envs/nextflow/' examples/sl4-vs-la2093/
# 期望：无输出（改为 nextflow26）
# 2026-08-06 注：本项在 2026-08-05 修复轮被遗漏，已于 2026-08-06 补修
```
