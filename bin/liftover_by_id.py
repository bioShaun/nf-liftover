#!/usr/bin/env python3
"""
根据 ID 文件（chrom_pos 格式）进行 liftover，生成目标基因组的 BED/ID 文件。
"""

from __future__ import annotations

from inspect import cleandoc
from pathlib import Path
from typing import Annotated
import subprocess

import pandas as pd
import typer
from loguru import logger
from pyfaidx import Fasta

MODULE_HELP = cleandoc(
    """
    根据 ID 文件（chrom_pos 格式）进行 liftover，生成目标基因组的 BED/ID 文件。

    \b
    使用示例:
      python liftover/liftover_by_id.py \\
        snp.id \\
        ref2query.chain \\
        ref.fa \\
        query.fa \\
        outdir

    \b
      强制重新生成:
      python liftover/liftover_by_id.py \\
        snp.id ref2query.chain ref.fa query.fa outdir --force

    \b
    输出格式:
      - <probe_name>.id: 1 列，pos_id（chrom_pos）
      - <probe_name>.bed: 3 列（chrom, start, pos），按 query.fa.fai 排序
      - <probe_name>.pos.tsv: 3 列（chrom, pos, id），其中 id 为输入 ID，chrom/pos 为 liftover 后坐标
      - <probe_name>.snpcalling.bed: 3 列（chrom, start, end），slop+merge 后的区间
    """
)

app = typer.Typer(help=MODULE_HELP, no_args_is_help=True)


# ---------------------------------------------------------------------------
# Pure functions
# ---------------------------------------------------------------------------


def parse_id_to_chrom_pos(snp_id: str) -> tuple[str, int]:
    """解析 'chrom_pos' 格式 ID，返回 (chrom, pos)。"""
    parts = snp_id.rsplit("_", maxsplit=1)
    if len(parts) != 2:
        raise ValueError(f"无法解析 ID: {snp_id!r}，预期格式 chrom_pos")
    chrom, pos_str = parts
    return chrom, int(pos_str)


def fetch_ref_nucleotide(ref_fa: Fasta, chrom: str, pos: int) -> str:
    """从参考基因组获取指定位置的碱基。"""
    return str(ref_fa[chrom][pos - 1 : pos].seq) if chrom in ref_fa else "N"


def make_id_vcf(id_file: Path, ref_fa_path: Path, *, force: bool = False) -> Path:
    """根据 ID 文件生成最小 VCF，用于 transanno liftvcf 输入。"""
    id_vcf_file = Path(f"{id_file}.vcf")
    if id_vcf_file.exists() and not force:
        logger.info(f"VCF 文件已存在，跳过生成: {id_vcf_file}")
        return id_vcf_file

    ref_fasta = Fasta(str(ref_fa_path))
    with open(id_file) as id_inf, open(id_vcf_file, "w") as vcf_out:
        vcf_out.write("##fileformat=VCFv4.2\n")
        vcf_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for line in id_inf:
            each_id = line.strip()
            if not each_id:
                continue
            chrom, pos = parse_id_to_chrom_pos(each_id)
            ref_seq = fetch_ref_nucleotide(ref_fasta, chrom, pos)
            vcf_out.write(f"{chrom}\t{pos}\t{each_id}\t{ref_seq}\t.\t.\t.\t.\n")

    logger.info(f"已生成 VCF: {id_vcf_file}")
    return id_vcf_file


def run_transanno_liftvcf(
    vcf: Path,
    chain: Path,
    ref_fa: Path,
    query_fa: Path,
    output: Path,
    rejected: Path,
) -> None:
    """调用 transanno liftvcf 进行 liftover。失败时抛出 RuntimeError。"""
    cmd = [
        "transanno",
        "liftvcf",
        "--original-assembly",
        str(ref_fa),
        "--new-assembly",
        str(query_fa),
        "--chain",
        str(chain),
        "--vcf",
        str(vcf),
        "--output",
        str(output),
        "--fail",
        str(rejected),
    ]
    logger.info(f"运行 liftover: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"transanno liftvcf 失败 (returncode={result.returncode})\n"
            f"stderr: {result.stderr}"
        )


def read_liftover_result(vcf_gz: Path) -> pd.DataFrame:
    """读取 liftover VCF 输出，返回去重的 DataFrame (chrom, pos, id, start, pos_id)。"""
    compression = "gzip" if vcf_gz.suffix == ".gz" else "infer"
    df = pd.read_table(
        vcf_gz,
        header=None,
        sep="\t",
        comment="#",
        usecols=[0, 1, 2],
        names=["chrom", "pos", "id"],
        compression=compression,
    )
    df["start"] = df["pos"] - 1
    df["pos_id"] = df["chrom"].astype(str) + "_" + df["pos"].astype(str)
    df.drop_duplicates(inplace=True)
    return df


def load_chrom_sizes_from_fai(fai: Path) -> pd.DataFrame:
    """从 .fai 文件加载染色体大小，返回 DataFrame (chrom, chrom_size)。"""
    chrom_df = pd.read_table(
        fai, header=None, names=["chrom", "chrom_size"], usecols=[0, 1]
    )
    chrom_df["chrom"] = chrom_df["chrom"].astype(str)
    chrom_df["chrom_size"] = pd.to_numeric(
        chrom_df["chrom_size"], errors="raise"
    ).astype(int)
    return chrom_df


def sort_bed_by_fai(bed_df: pd.DataFrame, chrom_df: pd.DataFrame) -> pd.DataFrame:
    """按 FAI 染色体顺序排序 BED DataFrame。"""
    df = bed_df.copy()
    df["chrom"] = pd.Categorical(
        df["chrom"].astype(str),
        categories=chrom_df["chrom"].tolist(),
        ordered=True,
    )
    return df.sort_values(by=["chrom", "start"]).reset_index(drop=True)


def slop_and_merge(
    bed_df: pd.DataFrame,
    chrom_sizes: dict[str, int],
    flank: int,
) -> pd.DataFrame:
    """
    对每个区间两侧扩展 flank bp（clamp 到 [0, chrom_size]），然后合并重叠区间。

    输入需要 chrom, start 列，以及 pos 或 end 列。
    返回 DataFrame (chrom, start, end)。
    """
    df = bed_df.copy()
    end_col = "end" if "end" in df.columns else "pos"
    df["end"] = df[end_col]

    df["start"] = (df["start"] - flank).clip(lower=0)
    df["end"] = df.apply(
        lambda row: min(row["end"] + flank, chrom_sizes.get(str(row["chrom"]), row["end"] + flank)),
        axis=1,
    )

    merged_rows: list[dict[str, object]] = []
    for chrom, group in df.groupby("chrom", sort=False, observed=True):
        intervals = group.sort_values("start")[["start", "end"]].values.tolist()
        if not intervals:
            continue
        cur_start, cur_end = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_end:
                cur_end = max(cur_end, e)
            else:
                merged_rows.append({"chrom": chrom, "start": int(cur_start), "end": int(cur_end)})
                cur_start, cur_end = s, e
        merged_rows.append({"chrom": chrom, "start": int(cur_start), "end": int(cur_end)})

    return pd.DataFrame(merged_rows, columns=["chrom", "start", "end"])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@app.command(help=MODULE_HELP)
def main(
    id_file: Annotated[Path, typer.Argument(help="输入 ID 文件（每行一个 chrom_pos 格式 ID）")],
    chain: Annotated[Path, typer.Argument(help="chain 文件路径")],
    ref_fa: Annotated[Path, typer.Argument(help="参考基因组 FASTA 路径")],
    query_fa: Annotated[Path, typer.Argument(help="目标基因组 FASTA 路径")],
    outdir: Annotated[Path, typer.Argument(help="输出目录")],
    force: Annotated[bool, typer.Option(help="强制重新生成所有中间文件")] = False,
    flank: Annotated[int, typer.Option(help="snpcalling BED 区间两侧扩展长度")] = 100,
) -> None:
    """根据 ID 文件进行 liftover 并生成 BED/ID 文件。"""
    outdir.mkdir(exist_ok=True, parents=True)

    # 1. 生成 VCF
    vcf = make_id_vcf(id_file, ref_fa, force=force)

    # 2. 运行 transanno liftvcf
    probe_name = id_file.stem
    lift_over_vcf = outdir / f"liftover.{id_file.name}.vcf.gz"
    rejected_vcf = outdir / f"rejected.{id_file.name}.vcf.gz"

    if force or not lift_over_vcf.is_file():
        run_transanno_liftvcf(vcf, chain, ref_fa, query_fa, lift_over_vcf, rejected_vcf)
    else:
        logger.info(f"Lifted VCF 已存在，跳过 liftover: {lift_over_vcf}")

    # 3. 读取结果
    lift_bed = read_liftover_result(lift_over_vcf)

    # 4. 按 FAI 排序
    query_fai = Path(f"{query_fa}.fai")
    chrom_df = load_chrom_sizes_from_fai(query_fai)
    sorted_bed = sort_bed_by_fai(lift_bed, chrom_df)

    # 5. 写入 probe_name.id
    probe_id_file = outdir / f"{probe_name}.id"
    sorted_bed.to_csv(
        probe_id_file, sep="\t", index=False, header=False, columns=["pos_id"]
    )
    logger.info(f"已写入 ID 文件: {probe_id_file}")

    # 6. 写入 probe_name.bed
    probe_bed = outdir / f"{probe_name}.bed"
    sorted_bed.to_csv(
        probe_bed, sep="\t", index=False, header=False, columns=["chrom", "start", "pos"]
    )
    logger.info(f"已写入 BED 文件: {probe_bed}")

    # 7. 写入 probe_name.pos.tsv
    probe_pos_file = outdir / f"{probe_name}.pos.tsv"
    sorted_bed.to_csv(
        probe_pos_file, sep="\t", index=False, header=False, columns=["chrom", "pos", "id"]
    )
    logger.info(f"已写入位点坐标文件: {probe_pos_file}")

    # 8. 生成 snpcalling BED (slop + merge)
    chrom_sizes_dict = dict(zip(chrom_df["chrom"], chrom_df["chrom_size"]))
    snpcalling_df = slop_and_merge(sorted_bed, chrom_sizes_dict, flank)
    snpcalling_sorted = sort_bed_by_fai(snpcalling_df, chrom_df)
    snpcalling_bed = outdir / f"{probe_name}.snpcalling.bed"
    snpcalling_sorted.to_csv(
        snpcalling_bed, sep="\t", index=False, header=False, columns=["chrom", "start", "end"]
    )
    logger.info(f"已写入 snpcalling BED: {snpcalling_bed}")


if __name__ == "__main__":
    app()
