# mm2-plus 替代 minimap2 加速基准测试 Spec

## 背景

[mm2-plus](https://github.com/at-cg/mm2-plus)（v1.3）是 minimap2 `2.31-r1302` 的加速分支，声称是 minimap2 的 drop-in 替代，输出近乎一致。其优化点：

1. 并行 chaining
2. 基于区间树的主链选择算法
3. AVX2/AVX512 SIMD 比对（来自 mm2-fast）
4. AVX2/AVX512 SIMD chaining（来自 Intel TAL）
5. 并行排序

论文报告基因组间比对最高 7x 加速（48 线程、完整基因组）。

nf-liftover 的 `ALIGN_WHOLE_CHROMOSOME` 是"整条染色体作为单条 query"的比对形态。minimap2 在这种形态下，chaining 阶段基本无法用满多线程（线程并行粒度是 query 序列），因此这正是 mm2-plus 并行 chaining 的目标场景。

## 目标

量化同参数下 minimap2 与 mm2-plus 的 wall time、CPU 利用率、内存峰值，并确认 PAF 输出一致性，产出"是否接入流程"的判定依据。

## 非目标

- 不修改 `modules/local/*.nf` 与 `conf/*.config`
- 不新增 `--aligner` 开关（仅约束本轮实验范围；已被后续小麦切窗口实验 supersede，最终决策见 [`docs/wheat-split-benchmark-spec.md`](wheat-split-benchmark-spec.md) 的最终结论）
- 不改动现有 `nf-liftover-tools` conda 环境
- 所有产物只落在 `/project/tmp/mm2plus-bench/`，不写入仓库

## 环境事实（已勘察）

| 项 | 值 |
| --- | --- |
| minimap2 | `/project/software/miniforge3/envs/nf-liftover-tools/bin/minimap2`，`2.31-r1302` |
| mm2plus | `/project/software/miniforge3/envs/mm2plus/bin/mm2plus`，`1.3`（bioconda `mm2plus=1.3`，已创建） |
| SIMD 分发 | 自动选 `mm2plus.avx512`（本机具备 avx512f/bw/dq/vl/vbmi/ifma） |
| stdout 纯净性 | 启动器提示（`Launching executable ...`）只写 stderr，stdout 干净，`> out.paf` 重定向安全 |
| 机器 | 32 核 / 122 GB RAM / gcc 15.2 / `/project/tmp` 所在盘剩余约 542 GB |
| 基准数据 | `IRGSP-1.0` 染色体 `1`：43,270,923 bp；`MH63RS3` 染色体 `Chr01`：45,027,022 bp |

关键前提：本机 minimap2 版本与 mm2-plus 的上游基线版本完全一致（均为 2.31-r1302），因此 PAF 可以逐行 diff，差异不会来自版本代差。

本机没有 `examples/sl4-vs-la2093/run.sh` 引用的番茄真实基因组（`/public/data/genomes/solanum_lycopersicum_LA2093/` 不存在），仓库自带 smoke 数据只有 10 kb，无法测速，故改用本地已有的水稻基因组。

## 数据准备

```bash
mkdir -p /project/tmp/mm2plus-bench && cd /project/tmp/mm2plus-bench
S=/project/software/miniforge3/envs/nf-liftover-tools/bin/samtools
$S faidx /public/data/genome/oryza_sativa_IRGSP-1.0/genome.fa 1     > irgsp.chr1.fa   # ref 侧
$S faidx /public/data/genome/oryza_sativa_MH63RS3/genome.fa Chr01   > mh63.chr1.fa    # query 侧
```

这两个 FASTA 已经提取完成，位于 `/project/tmp/mm2plus-bench/`（`irgsp.chr1.fa` 42 MB、`mh63.chr1.fa` 44 MB），可直接跑基准。所有中间产物统一放 `/project/tmp`，不使用 `/tmp`。

参照 `modules/local/align_whole_chromosome.nf`：

```
minimap2 ${args} -t ${task.cpus} "${pair_fastas.query_chrom_fa}" "${pair_fastas.ref_chrom_fa}"
```

即 **target = query 基因组染色体（mh63），query = ref 基因组染色体（irgsp）**。顺序不可颠倒，否则与流程行为不一致，耗时和输出都不可比。

## 基准方法

严格复用 `conf/modules.config` 中 `ALIGN_WHOLE_CHROMOSOME` 的 `ext.args = '-cx asm5 --cs'`。线程档位取两个：

- `16`：`conf/base.config` 的 `large_mem`（`params.max_cpus = 16`），即 `ALIGN_WHOLE_CHROMOSOME` 实际使用值
- `8`：`split_mem`（`ALIGN_SPLIT_WINDOW` / `BUILD_QUERY_MMI` 使用值）

两个工具**串行**执行，避免互相抢占 CPU 影响计时。

`/project/tmp/mm2plus-bench/bench.sh`：

```bash
#!/usr/bin/env bash
set -euo pipefail

cd /project/tmp/mm2plus-bench

MM2=/project/software/miniforge3/envs/nf-liftover-tools/bin/minimap2
MM2P=/project/software/miniforge3/envs/mm2plus/bin/mm2plus
ARGS='-cx asm5 --cs'
TARGET=mh63.chr1.fa
QUERY=irgsp.chr1.fa

for t in 16 8; do
  for tool in mm2 mm2plus; do
    bin=$MM2
    [ "$tool" = mm2plus ] && bin=$MM2P
    out="${tool}.t${t}.paf"
    log="${tool}.t${t}.time"
    echo "=== running ${tool} -t ${t} ==="
    /usr/bin/time -v "$bin" $ARGS -t "$t" "$TARGET" "$QUERY" > "$out" 2> "$log" || {
      echo "FAILED ${tool} t=${t}"
      tail -20 "$log"
      exit 1
    }
    grep -E 'Elapsed \(wall clock\)|Percent of CPU|Maximum resident' "$log"
    wc -l "$out"
  done
done

echo "=== done ==="
```

执行顺序为 `mm2 t16 → mm2plus t16 → mm2 t8 → mm2plus t8`。单次预计数分钟量级，建议后台跑并轮询：

```bash
cd /project/tmp/mm2plus-bench && chmod +x bench.sh
nohup bash bench.sh > bench.log 2>&1 &
tail -f bench.log
```

需要回传的数据：每个 `*.time` 中的 `Elapsed (wall clock)`、`Percent of CPU this job got`、`Maximum resident set size`，以及每个 PAF 的行数。

## 正确性校验

```bash
cd /project/tmp/mm2plus-bench

# 1) 同线程档位下两工具是否完全一致
diff -q mm2.t16.paf mm2plus.t16.paf
diff -q mm2.t8.paf  mm2plus.t8.paf

# 2) 若不一致，统计差异规模
diff mm2.t16.paf mm2plus.t16.paf | grep -c '^[<>]'

# 3) 汇总指标对比：记录数、比对块总长、整体 identity
for f in mm2.t16.paf mm2plus.t16.paf; do
  awk -v f="$f" '{m+=$10; b+=$11} END{printf "%s\trecords=%d\taln_len=%d\tidentity=%.6f\n", f, NR, b, m/b}' "$f"
done

# 4) 线程非确定性对照：同工具不同线程数应输出一致
diff -q mm2.t16.paf    mm2.t8.paf
diff -q mm2plus.t16.paf mm2plus.t8.paf
```

第 4 步用于区分"工具间差异"与"线程调度导致的非确定性"。若同工具换线程数就已经不一致，则工具间的 diff 不能直接归因于 mm2-plus。

## 结果记录模板

| 工具 | 线程 | wall time | CPU% | Max RSS | PAF 行数 | 与 minimap2 diff 行数 |
| --- | --- | --- | --- | --- | --- | --- |
| minimap2 2.31 | 16 | | | | | — |
| mm2plus 1.3 | 16 | | | | | |
| minimap2 2.31 | 8 | | | | | — |
| mm2plus 1.3 | 8 | | | | | |

加速比 = `minimap2 wall time / mm2plus wall time`（分别按线程档位计算）。

汇总指标：

| 工具 | records | aln_len | identity |
| --- | --- | --- | --- |
| minimap2 2.31 | | | |
| mm2plus 1.3 | | | |

## 判定标准

- **通过（可考虑接入）**：`t16` 下 PAF diff 为空，或差异行占比 < 0.1%；wall time 加速 ≥ 1.3x；Max RSS 不超过 `large_mem` 上限（60 GB）。
- **需追加验证**：加速明显但 PAF 存在差异 → 追加水稻 12 对染色体全量比对，并在 chain 层面（`transanno minimap2chain`）对比结果差异是否影响 liftover 坐标。
- **不接入**：加速 < 1.1x，或输出差异会改变 chain 结构。

## 风险与备注

- 单条 43 Mb 染色体不足以复现论文级 7x（论文用 48 线程 + 完整基因组），`t16` 的结果偏保守，应作为下界看待。
- IRGSP-1.0（japonica）与 MH63RS3（indica）分化约 1%，`asm5` 偏严格；但只有与流程默认参数保持一致，测速结果才对本流程有意义。
- mm2-plus 官方说明：大基因组的 PAF 可能存在可忽略的细微差异；读长比对若要严格一致需加 `--max-chain-skip=1000000`。本基准是基因组间比对，不加该参数，保持与流程一致。
- bioconda 的 `mm2plus` 通过启动器按 CPU 特性分发到 `mm2plus.avx512`，无需自行编译；如后续发现分发失败，可直接调用 `mm2plus.avx512`。

## 若后续决定接入流程（仅备注，本次不实施）

改动面很小：

- `modules/local/build_query_mmi.nf`、`align_whole_chromosome.nf`、`align_split_window.nf` 三处的二进制名与 `versions.yml` 采集
- `nextflow.config.example` 的 `params.conda_envs` 增加一个 aligner 环境项，`conf/conda.config` 增加对应 label
- `conf/modules.config` 的 `ext.args` 可原样复用（mm2-plus 命令行与 minimap2 完全兼容）

## 执行结果（2026-07-27）

原始产物在 `/project/tmp/mm2plus-bench/`（`bench.sh` / `bench2.sh` / `*.time` / `*.paf` / `RESULTS.md`）。

### 主基准（`-cx asm5 --cs`，IRGSP chr1 → MH63 Chr01）

| 工具 | 线程 | wall time | CPU% | Max RSS | PAF 行数 | 加速比 |
| --- | --- | --- | --- | --- | --- | --- |
| minimap2 2.31 | 16 | 57.42s | 101% | 2.19 GiB | 1432 | — |
| mm2plus 1.3 | 16 | 47.75s | 121% | 2.33 GiB | 1432 | 1.20x |
| minimap2 2.31 | 8 | 57.77s | 100% | 2.20 GiB | 1432 | — |
| mm2plus 1.3 | 8 | 46.37s | 116% | 2.24 GiB | 1432 | 1.25x |

汇总指标两者完全相同：`records=1432`、`match=35369367`、`aln_len=41472591`、`identity=0.852837`。

### 正确性

- 同工具 `t16` vs `t8`：byte-identical，无线程非确定性。
- 跨工具同线程：原始 diff 有 72 处差异标记，但仅为记录顺序不同；`sort` 后 byte-identical。
- chain 层（`transanno minimap2chain`）：内容在忽略 chain id 编号与记录顺序后完全一致（1432 条全部对应），但原始文件 diff 552 行，因为 chain id 按输入顺序递增分配。

### 阶段分解（追加实验，t16）

| 阶段 | minimap2 | mm2plus | 加速比 |
| --- | --- | --- | --- |
| 仅 mapping（`-x asm5`，无 `-c`） | 13.26s / CPU 105% | 9.23s / CPU 227% | 1.44x |
| 完整（`-cx asm5 --cs`） | 57.42s / CPU 101% | 47.75s / CPU 121% | 1.20x |
| 推算 alignment 段 | ~44.2s | ~38.5s | ~1.15x |

两者 total CPU time 几乎相同（58.1s vs 57.9s），说明 mm2-plus 并未减少总计算量，1.20x 完全来自把 chaining 段摊到其他线程。chaining 段确实并行有效（CPU 105% → 227%，1.44x），但它只占整体约 23%；剩余约 77% 是单条 query 的碱基级比对，只有约 1.15x。因此**加线程无收益（t8 ≈ t16），该形态的加速上限本身就低**。

### split 模式对照（追加实验，同数据切 5 × 10 Mb 窗口，单进程 t16）

| 场景 | minimap2 | mm2plus |
| --- | --- | --- |
| 整条染色体 | 57.42s / CPU 101% | 47.75s / CPU 121% |
| 5 × 10 Mb 窗口 | 17.27s / CPU 355% | 13.47s / CPU 389% |

切窗口本身带来约 3.3x，远大于换 aligner 的 1.20x。注意窗口版输出 1473 条 vs 整条 1432 条（边界效应，需 `restore_split_paf.py` 还原），Max RSS 升至约 4.8 GB，并非等价替换。当前 `split_threshold = 100_000_000`，水稻（43 Mb）与番茄（约 90 Mb）染色体都走 whole 模式，故这是配置调优议题，与 mm2-plus 无关。

## 结论

**暂不接入 mm2-plus。**

| 检查项 | 结果 |
| --- | --- |
| PAF 一致性 | 通过（内容级完全一致） |
| Max RSS ≪ 60 GB | 通过（约 2.3 GiB） |
| wall time ≥ 1.3x | 不通过（1.20x） |

1.20x 落在规范的 1.1–1.3 中间带。阻碍因素不是"并行 chaining 没吃满"，而是运行时间的约 77% 花在单 query 碱基级比对上，mm2-plus 在该段只有约 1.15x，受 Amdahl 限制整体难以突破。

后续动作修正：

- **不建议**按原规范跑水稻 12 对染色体，同物种同分化度只会把 1.20x 重复 12 次，阶段占比不变，信息量低。
- 若仍要争取更高上界，应换分化度更高 / 染色体更大的真实对（例如番茄 SL4 vs LA2093 约 90 Mb 且分化更大），此时 chaining 段占比上升，mm2-plus 收益才可能变大。
- 若目标只是缩短流程 wall time，优先评估 `split_threshold` / `split_size` 调优（同数据 3.3x），而不是换 aligner。
- 若将来确定接入，需注意 chain 文件跨 aligner 不可 byte 比较（chain id 编号与记录顺序会变），任何对 `all.chain` 的 golden-file 或哈希断言必须先归一化。
