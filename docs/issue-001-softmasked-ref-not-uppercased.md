# Issue 001 — 软屏蔽参考下 `--id` 模式的 REF 未大写，transanno 全判 `UNEXPECTED_REF`

- **状态**：已修复（2026-07-30，`fetch_ref_nucleotide` 大写化 + 单测 `test_fetch_ref_nucleotide_uppercases_softmasked_bases`）
- **发现时间**：2026-07-30
- **发现场景**：TC-BY-peanut，Fuhuasheng (GCA_004170445.1) → Tifrunner gnm1 liftover
- **影响面**：所有用 **softmasked 参考** 跑 `--id` 模式的项目
- **严重度**：高 —— 静默丢位点，日志里只表现为"一部分位点 liftover 失败"，看不出是 bug
- **修复成本**：一行

## 现象

`/data_0/liftover/arachis_hypogaea_Fuhuasheng-to-arachis_hypogaea_tifrunner/TC-BY-Peanut-Fuhuasheng/`

| | 位点数 |
|---|---|
| 输入 `.id` | 130,328 |
| `liftover.*.id.vcf.gz` 成功 | 105,365 |
| `rejected.*.id.vcf.gz` 失败 | 24,963（`FAILED_REASON=UNEXPECTED_REF`） |

按 REF 列的大小写拆开看，规律非常干净：

| | REF 小写 | REF 大写 |
|---|---|---|
| 成功 105,365 条 | **0** | 105,365 |
| 失败 24,963 条 | **18,020** | 6,943 |

同一批流程跑的 NDH108 → Tifrunner（参考 `NDH108.fa` **不是** softmasked）1,655 条失败**全是大写**，
不受影响 —— 这条对照排除了"重复区本来就比不上"的解释。

## 根因

`bin/liftover_by_id.py:63-65`：

```python
def fetch_ref_nucleotide(ref_fa: Fasta, chrom: str, pos: int) -> str:
    """从参考基因组获取指定位置的碱基。"""
    return str(ref_fa[chrom][pos - 1 : pos].seq) if chrom in ref_fa else "N"
```

`pyfaidx` 原样返回 FASTA 里的字符。参考若是 softmasked，重复区就是小写 `a/c/g/t`，
`make_id_vcf()` 把它直接写进 VCF 的 REF 列（`liftover_by_id.py:82`）。
`transanno liftvcf` 校验 REF 时按**字面**比对，`c` != `C`，于是整条判 `UNEXPECTED_REF`。

顺带一个次要问题：同一行的 `else "N"` 兜底（contig 不在参考里）也会被 transanno 判
`UNEXPECTED_REF`，跟"碱基对不上"混在同一个失败原因里，排查时分不开。

## 复现 / 验证实验

从 `rejected.*.vcf.gz` 里抽 2,000 条 REF 为小写的记录，做成两份**除大小写外逐字节相同**的 VCF，
用同一条 `all.chain`、同一个 `-r/-q` 参考跑 `transanno liftvcf`：

| VCF | 通过 | 失败 |
|---|---|---|
| 原样（小写 REF） | **0** | 2,000 |
| 只把 REF 转大写 | **1,870（93.5%）** | 130 |

脚本（可直接拿来复跑）：
`/data_0/panel_design/projects/TC-BY-peanut/design/runs/DR-20260730-001/config/make_lowercase_ref_probe.py`
`/data_0/panel_design/projects/TC-BY-peanut/design/runs/DR-20260730-001/config/run_chain_liftover_crosscheck.sh`

按 93.5% 外推，Fuhuasheng 这一轮**白丢约 16,850 个位点**（18,020 × 0.935），
占输入的 13%，占全部失败的 68%。剩下 6,943 条大写失败的才是真的比不上。

## 建议修复

`bin/liftover_by_id.py`：

```python
def fetch_ref_nucleotide(ref_fa: Fasta, chrom: str, pos: int) -> str:
    """从参考基因组获取指定位置的碱基（softmasked 参考会给小写，必须大写后再写进 VCF，
    否则 transanno 按字面比对 REF 会全判 UNEXPECTED_REF）。"""
    if chrom not in ref_fa:
        return "N"
    return str(ref_fa[chrom][pos - 1 : pos].seq).upper()
```

配套建议（都不是必须，但能防止同类问题再静默发生）：

1. **建完 VCF 就自查**：统计 REF 列里非 `ACGTN` 的比例，>0 直接 `logger.warning` 并把
   softmasked 参考这件事说出来。
2. **`N` 兜底与真失败分开**：contig 不在参考里的位点应该在建 VCF 阶段就单独列出来，
   不要塞进 VCF 让 transanno 报 `UNEXPECTED_REF`。
   → 已拆为 issue-002（`docs/issue-002-missing-contig-sites-not-separated.md`）跟进。
3. **收尾报告加通过率断言**：`--id` 模式下通过率低于某个阈值（比如 90%）时显式告警。
   这一轮 81% 的通过率在日志里没有任何提示。

## 备注

- 该项目的 Fuhuasheng 目录下有 `realign_blast/` 阶段的日志，但最终 `.pos.tsv` 的行数
  （105,365）与 transanno 成功数完全相等，说明这 ~1.8 万条并没有被 blast rescue 捞回来。
  即使能捞回来，用 BLAST 重比对去补一个 `.upper()` 就能解决的问题，代价也差了一个量级。
- 完整的复核过程（含另一条与本 issue 无关的、关于 BLAST 侧翼投影加区间约束的结论）写在
  `/data_0/panel_design/projects/TC-BY-peanut/docs/plan-2026-07-30-functional-gene-and-customer-sites.md` §A6'。
