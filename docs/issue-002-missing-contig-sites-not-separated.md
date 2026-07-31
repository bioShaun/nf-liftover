# Issue 002 — contig 不在参考里的位点应单独列出，不要塞进 VCF 混进 `UNEXPECTED_REF`

- **状态**：已修复（2026-07-30，`make_id_vcf()` 分流缺失 contig 位点 + 单测
  `test_make_id_vcf_separates_sites_with_missing_contigs` /
  `test_make_id_vcf_raises_when_all_contigs_missing`）
- **来源**：issue-001 配套建议第 2 条（2026-07-30 拆出跟进）
- **影响面**：所有 `--id` 模式的项目；位点 ID 的 contig 名与参考 FASTA 不匹配时触发
- **严重度**：中 —— 不丢位点（反正也 liftover 不了），但失败原因被污染，排查时分不开
- **修复成本**：小 —— 集中在 `bin/liftover_by_id.py` 的 `make_id_vcf()`

## 现象

`bin/liftover_by_id.py` 的 `fetch_ref_nucleotide()` 对 `chrom not in ref_fa` 的位点兜底返回 `"N"`，
`make_id_vcf()` 照常在 VCF 里写一行 `REF=N`。transanno 拿 `N` 去和参考比对，必然判
`UNEXPECTED_REF`——和"碱基真的对不上"（如 issue-001 的 softmasked 问题）混在同一个
`FAILED_REASON` 里。

后果：

- rejected VCF 里 `UNEXPECTED_REF` 的数量不能用来评估"参考/位点不匹配程度"，
  必须先把 contig 缺失的部分剔掉才知道真失败有多少。
- contig 名写错（比如 ID 用 `chr1` 而参考用 `SL4.0ch01`）这种**全量性输入错误**，
  只表现为"所有位点 liftover 失败"，日志里没有任何针对性的提示。

## 建议修复

在 `make_id_vcf()` 建 VCF 阶段就把两类位点分开：

1. `chrom not in ref_fa` 的位点**不写进 VCF**，单独输出到
   `<id_file>.missing_contigs.tsv`（或类似命名），列出 id、chrom、pos。
2. 有缺失时 `logger.warning` 明确报数量，并提示检查 ID 的 contig 命名是否与参考一致。
3. 全部为缺失（即大概率 contig 命名整体不匹配）时考虑直接报错退出，
   而不是让 transanno 空跑一遍再给一个 100% 失败的 rejected VCF。

这样 transanno 输出的 `UNEXPECTED_REF` 就只剩"位点在参考里但碱基对不上"的真实失败。

## 备注

- issue-001 已修掉 softmasked 小写 REF 的问题；issue-001 修复后 `UNEXPECTED_REF`
  的剩余构成 = 真比对失败 + contig 缺失兜底（本 issue）。
- 配套建议第 1 条（建完 VCF 统计 REF 列非 `ACGTN` 比例并告警）和第 3 条
  （通过率低于阈值告警）未拆 issue，如需要可再补。
