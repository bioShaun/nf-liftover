# Issue 003 — BLAST rescue 的验收闸门抓不住真错，且比它要"救"的 transanno 更宽松

- **状态**：已修复（2026-07-30，tc-pytools `realign_blast.py` 与 nf-liftover `merge_blast_rescue.py`
  两侧均已落地，详见文末"修复记录"）
- **发现时间**：2026-07-30
- **发现场景**：TC-BY-peanut，用 BLAST 侧翼投影做 Fuhuasheng / NDH108 → Tifrunner liftover
  后，与 `/data_0/liftover` 的 minimap2 → chain → transanno 结果全量交叉复核（2034 个位点）
- **影响面**：所有走 `merge_blast_rescue.py` 捞回 transanno 失败位点的项目；
  多倍体 / 有近期全基因组复制的物种尤其严重
- **严重度**：高 —— 静默输出错误坐标。四个 reject code 全过，但坐标落在旁系同源拷贝上
- **修复成本**：中 —— 第 1、3、4 条需要动 `tc-pytools/panel/realign_blast.py`（上游产 `selection.tsv` 的地方）

## 代码链条

```
tc-pytools/panel/realign_blast.py   ← 真正跑 blastn、选 hit、推坐标
        ↓  *.selection.tsv
nf-liftover/bin/merge_blast_rescue.py  ← 只做验收和合并
```

附带说明：`merge_blast_rescue.py` 目前**没有被任何 `.nf` 引用**
（`grep -rl blast modules/local subworkflows conf main.nf nextflow.config` 无匹配），
是手工脚本。TC-BY-peanut 那一轮的 `.pos.tsv` 行数（105,365）与 transanno 成功数完全相等，
说明 rescue 结果从未被并入交付。要正式接进流水线，建议先解决下面四条。

## 实测背景数据

用一条独立的 `all.chain` 重跑 transanno 覆盖全部 2034 个候选位点，与 BLAST 投影逐位点比对，
得到 284 个分歧位点。对每个分歧位点**去掉染色体/区间约束重跑 BLAST** 做仲裁：

| 分歧类别 | n | 受约束命中 pident 中位 | 自由最优命中 pident 中位 | 自由命中站 chain 一边 |
|---|---|---|---|---|
| `disagree_chrom`（跨染色体，**真错**） | 113 | **87.82** | 100.00 | **113 / 113** |
| `disagree_pos`（同染色体位置差，**真错**） | 56 | 94.05 | 100.00 | 50 / 56 |
| `disagree_homoeolog`（A/B 拷贝，**真歧义**） | 115 | 100.00 | 100.00 | 15 / 115 |

下面每条问题都拿这批数字量化。

---

## 问题 1（最高优先）：`match_ratio` 是覆盖度不是一致度，`pident` 被丢掉了

`tc-pytools/panel/realign_blast.py:775`：

```python
df["match_ratio"] = df["align_len"] / df["informative_len"]
```

`align_len` 是 BLAST 的比对长度，**错配不会让它变短**。一条 88% 一致度、全长比上的命中，
`match_ratio ≈ 1.0`，轻松穿过 0.95。

更关键的是：`pident` 在 `LEGACY_BLAST_COLUMNS`（`:36`）里解析出来了，但
`ALIGNMENT_COLUMNS`（`:51-67`）没带它，`parse_blast_results` 返回时 `[ALIGNMENT_COLUMNS]` 把它切掉。
于是**一致度既进不了筛选，也进不了 `selection.tsv`**
（`tests/test_core_helpers.py:328` 的表头可证：有 `match_ratio` / `align_len` / `mismatches`，无 `pident`）。

`merge_blast_rescue.py` 的 `--min-match-ratio 0.95` 完整继承了这个错觉：
参数名像"95% 匹配"，实际含义是"95% 覆盖"，对旁系同源假命中几乎没有区分力。

**实测**（TC-BY-peanut 仲裁轮的原始 blastn 输出，1018 条命中）：

| pident 在 84–95% 之间的命中 | n = 346 |
|---|---|
| query 覆盖度（`length / qlen`）中位数 | **0.978** |
| 覆盖度 > 0.95 | 211 条（**61.0%**） |
| 覆盖度 > 0.90 | 231 条（66.8%） |

**六成的旁系同源假命中直接穿过 0.95 闸门。**
反过来，若存在一条 95% 的一致度底线，那 113 个 `disagree_chrom` 错误里 **97 个会被拦下**
（受约束命中只有 16/113 达到 pident ≥ 95）。

### 建议

1. `ALIGNMENT_COLUMNS` 加回 `pident`，`SELECTION_REPORT_COLUMNS` 输出它 ——
   这是判断旁系同源误置最直接的一个数（87.8 vs 100.0），交付表里现在完全看不到。
2. `build_id_mapping()` 加 `min_pident` 参数，默认给一个保守值（多倍体建议 ≥ 95）。
3. `merge_blast_rescue.py` 增加 `--min-pident` 与对应 reject code `BLAST_LOW_IDENTITY`；
   同时考虑把 `--min-match-ratio` 改名为 `--min-query-coverage`（保留旧名做 alias），
   避免继续误导。
   - 过渡期内可用 `1 - mismatches / align_len` 近似，`selection.tsv` 里这两列都有。

## 问题 2：`--id-chrom-map` 的硬过滤会把最优命中删掉且不留痕迹

`realign_blast.py:674-705` 的 `_filter_by_id_chrom_map()` 对**登记过的 id 硬删**
非允许染色体上的 hit。关键在顺序：

```python
df = _filter_by_id_chrom_map(df, id_chrom_map)   # :782  先删
...
df["rank"] = df.groupby("id", sort=False).cumcount() + 1   # :818  后排名
```

`rank` 和 `selection_reason` 都是在**被截断过的候选集**上生成的。被删掉的那条 100% 命中
在任何输出里都不留痕迹 —— `selection.tsv` 只会写
`rank #1; bitscore=…; id_chrom_map: 严格限制通过`，读起来像是干净的最优命中。

日志侧只有"未登记 id"的 warning（`:693`）。**"登记了、但把最优命中删掉了"这个真正危险的
情况一声不吭。**

这与本项目 BLAST 投影里 `--expect-interval` 的失效模式完全同构：一旦上游把位点桥接到了
**旁系同源**基因，约束就会强行选中一个低一致度命中，把基因层面的错误传导成坐标错误。
实测代价见上表：`disagree_chrom` 113 个位点，自由最优命中 **113/113 站独立 chain 一边、
0/113 站受约束结果一边**。

`--chr-map`（`:795-808`）温和一些 —— 它是排序键而非过滤器，无匹配时会回退到跨染色体最优。
但仍然是"匹配但更差"压过"不匹配但更好"，而且 `chr_map_status == "matched"`
分不出"匹配且本来就是全局最优"和"只因为被强推才排第一"。
`merge_blast_rescue.py` 的 `BLAST_CHR_MAP_NOT_MATCHED` 建立在这个分不开的状态上。

### 建议

1. **过滤前先记录全局最优**：在 `_filter_by_id_chrom_map()` 之前算出每个 id 的
   `global_best_bitscore` / `global_best_chrom` / `global_best_pident`，作为列带到
   `selection.tsv`。审查时一眼能看出"被约束换掉了什么"。
2. 新增状态值区分：`strict_and_global_best` / `strict_but_dropped_better_hit`；
   `chr_map_status` 同理加 `matched_but_not_global_best`。
3. `merge_blast_rescue.py` 增加 reject code `BLAST_CONSTRAINT_OVERRODE_BEST_HIT`：
   被约束选中的命中 pident 比全局最优低超过阈值（建议 2 个百分点）时不予采信。
4. 汇总日志里报出"因约束而更换最优命中"的位点数 —— 这个数偏高就说明上游的
   id→染色体映射本身有问题，值得先回头查桥接而不是继续 rescue。

## 问题 3：`_has_tied_top` 用严格相等，正好挑错了要抓的那一类

`bin/merge_blast_rescue.py` 的歧义闸门：

```python
top_bitscore = float(top.iloc[0]["bitscore"])
return (group["bitscore"].astype(float) == top_bitscore).sum() > 1
```

浮点等值。对照实测的 best-vs-second bitscore 差：

| 类别 | n | margin 中位 | 严格并列（margin = 0） | ≤ 15 bits |
|---|---|---|---|---|
| `disagree_chrom`（**真错，可修**） | 113 | 11.0 | **0 / 113** | 58 / 113 |
| `disagree_pos`（**真错，可修**） | 56 | 50.5 | 1 / 56 | 3 / 56 |
| `disagree_homoeolog`（**真歧义，不可修**） | 115 | 0.0 | 69 / 115 | 102 / 115 |

**这个闸门一个真错都抓不到（0/113），却把不可解的 A/B 拷贝歧义拦下了大半。**
它精确地做反了：该放行的拦了，该拦的放了。

此外还有一个结构性盲区：`_has_tied_top` 读的 `selection.tsv` 已经是
**id_chrom_map 过滤后 + top-N 截断后**的结果。如果问题 2 的硬过滤删掉了真最优，
margin 检查根本看不见对手，必然判"无并列"。

### 建议

1. 换成 margin 阈值而非等值：`best - second < max(20, 0.05 * best)` 即标记
   `BLAST_AMBIGUOUS_TOP_HIT`。阈值应可配。
2. 歧义判定必须在**任何约束过滤之前**做，或至少把过滤前的
   `global_second_best_bitscore` 带进 `selection.tsv` 供其使用（与问题 2 的建议 1 共用）。
3. margin ≈ 0 的位点在多倍体里是**真歧义、两种方法都无法证伪**，
   正确处理是标注后交给下游/客户判断，不要假装解决了 ——
   本项目最终就是把这 115 个位点单独标成 `homoeolog_ambiguous` 交付的。
4. `--report-top-n` 默认 3 时 `selection.tsv` 才有次优可比；
   `merge_blast_rescue.py` 应在 group 内只有 1 条记录时显式报"无法判断歧义"，
   而不是当作"无并列"通过（当前 `_has_tied_top` 在这种情况下返回 `False`，即放行）。

## 问题 4：rescue 路径不核对 REF 碱基，比它要"救"的 transanno 还宽松

transanno 是**校验 REF 的**（issue-001 正是因此暴露）。而：

- `realign_blast.py` 算出 `pos`（`_infer_position`，`:429-435`）就结束，不看目标基因组上
  那个碱基是什么；
- `merge_blast_rescue.classify_blast_rescue()` 的四个检查
  （`BLAST_LOW_MATCH_RATIO` / `BLAST_CHR_MAP_NOT_MATCHED` /
  `BLAST_TARGET_CHROM_NOT_IN_MAPPING` / `BLAST_TIED_TOP_HIT`）也没有一个碰目标碱基。

结果：**这条 rescue 可以把 transanno 正确拒绝的位点重新放进交付**，
而且是在"比原方法更弱的证据"下放进去的。

本项目的 BLAST 投影是走比对串定位到中心碱基后，**核对 Tifrunner 实际碱基与 REF/ALT**
才输出坐标的。这一步很便宜（`samtools faidx` 批量取区间，2000 个区间一次调用），
相对 BLAST 本身的开销可忽略。

### 建议

1. `merge_blast_rescue.py` 增加必需参数 `--target-fasta`，对每条待采信的 rescue
   取出目标碱基，与源位点的 REF/ALT 比对，不符则 `BLAST_REF_MISMATCH`。
   （允许 REF↔ALT 互换的情况单独标 `BLAST_REF_IS_ALT`，这是真实且有意义的一类。）
2. 顺带能防住软屏蔽问题复发：取到小写碱基时先 `.upper()` 再比（见 issue-001）。
3. 汇总报告里给出 rescue 的 REF 核对通过率 —— 这个数低于 95% 就说明整批 rescue
   不可信，应该整体丢弃而不是逐条挑。

---

## 相关

- issue-001 / issue-002：`--id` 模式建 VCF 的两个问题（均已修复）
- `tc-pytools` 侧另有两条与本 issue 相邻、但属于那个仓库的问题
  （`.pos.tsv` / `.idmap.tsv` 因 `max_hits` 默认 3 变成一对多；
  `-max_target_seqs 10` + 默认 megablast 对异源四倍体偏紧），
  见 `tc-pytools/docs/realign_blast_output_contract_and_search_scope.md`
- 完整的交叉复核过程与数据：
  `/data_0/panel_design/projects/TC-BY-peanut/docs/plan-2026-07-30-functional-gene-and-customer-sites.md` §A6'
- 复核脚本（可直接复跑）：
  `/data_0/panel_design/projects/TC-BY-peanut/design/runs/DR-20260730-001/config/crosscheck_chain_liftover.py`

## 备注：写得对的部分

审查中逐项验算过、**没有问题**、改动时不要动的地方：

- BTOP 路径的坐标推算（`realign_blast.py:379-426`）与无 gap 路径（`:319-339`）在正负链两种
  约定下结果一致，无 off-by-one。BLAST tabular 里 query 坐标恒为升序、链向体现在 subject 上，
  两条路径分别用 `offset_fwd` 和 `offset_fwd/offset_rev + query_start` 重算，代数上等价。
- `informative_len = query_len - n_count`（`:759`）把 N 从 `match_ratio` 分母里扣掉，
  并对全 N 的 query 显式 warn 后丢弃 —— 比常见实现细致。
- `_iupac_to_representative()`（`:449-466`）保留 N 不展开成 A，避免 N→A 造成的假阳性匹配。


## 修复记录（2026-07-30）

**tc-pytools 侧**（`panel/realign_blast.py`，40 个单测全绿）：

- `pident` 加回 `ALIGNMENT_COLUMNS` / `SELECTION_REPORT_COLUMNS`，进入 `selection.tsv`。
- `build_id_mapping()` 新增 `min_pident=95.0`，CLI 两个命令暴露 `--min-pident`。
- 约束过滤前新增 `_annotate_global_best()`，`selection.tsv` 带
  `global_best_bitscore` / `global_best_chrom` / `global_best_pident` /
  `global_second_best_bitscore` 四列（口径：质量闸门之后、约束过滤之前）。
- `id_chrom_status`：`strict` 细分为 `strict_and_global_best` / `strict_but_dropped_better_hit`；
  `chr_map_status` 新增 `matched_but_not_global_best`；旧值 `fallback` / `n/a` 语义不变。
- 汇总日志报"因约束更换最优命中"的位点数。

**nf-liftover 侧**（`bin/merge_blast_rescue.py`，37 个单测全绿）：

- 新增 `--min-pident`（默认 95.0）→ `BLAST_LOW_IDENTITY`；无 `pident` 列时用
  `100 * (1 - mismatches / align_len)` 近似并 warning。
- `--min-match-ratio` 改名 `--min-query-coverage`（旧名保留为 alias；
  reject code `BLAST_LOW_MATCH_RATIO` 不变以免破坏旧报告解析）。
- 新增 `BLAST_CONSTRAINT_OVERRODE_BEST_HIT`：选中命中比 `global_best_pident` 低超过
  `--max-constraint-pident-drop`（默认 2.0）时拒收；旧格式无 `global_best_*` 列时降级跳过。
- `_has_tied_top` 浮点等值改为 margin 闸门：`best - second < max(--tied-margin=20,
  --tied-margin-frac=0.05 * best)` → `BLAST_AMBIGUOUS_TOP_HIT`（取代 `BLAST_TIED_TOP_HIT`）；
  优先用 `global_second_best_bitscore`；单记录且无法判断歧义时显式 warning。
- 新增必需参数 `--target-fasta` 做 REF 核对（pyfaidx 取碱基先 `.upper()`，防 issue-001 复发）：
  通过 / `BLAST_REF_IS_ALT`（单独统计）/ `BLAST_REF_MISMATCH`；汇总报告给
  `ref_check_pass_rate`。REF/ALT 来源：selection `alleles` 列 → `--rejected-vcf` 回退。
- `chr_map_status == "matched_but_not_global_best"` 视为"匹配"，是否采信交给
  pident 落差闸门判断（不再被 `BLAST_CHR_MAP_NOT_MATCHED` 一刀切）。

**仍未处理**（不在本 issue 范围）：`merge_blast_rescue.py` 尚未接进任何 `.nf` 流程
（接入时记得 `--target-fasta` 为必需参数）；tc-pytools 的 `max_hits` 默认 3 一对多、
`-max_target_seqs 10` 偏紧两个问题见 tc-pytools 自己的文档。
