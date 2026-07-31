#!/usr/bin/env python3
"""
根据 ID 文件（chrom_pos 格式）进行 liftover，生成目标基因组的 BED/ID 文件。
"""

from __future__ import annotations

from inspect import cleandoc
import subprocess
from pathlib import Path
from typing import Annotated, Sequence

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
      - 启用 split 时，.id/.pos.tsv 仍为 liftover 后原目标坐标，.bed/.snpcalling.bed 输出 split 坐标
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
    """从参考基因组获取指定位置的碱基（softmasked 参考会给小写，必须大写后再写进 VCF，
    否则 transanno 按字面比对 REF 会全判 UNEXPECTED_REF）。"""
    if chrom not in ref_fa:
        return "N"
    return str(ref_fa[chrom][pos - 1 : pos].seq).upper()


def make_id_vcf(id_file: Path, ref_fa_path: Path, missing_out: Path | None = None) -> Path:
    """根据 ID 文件生成最小 VCF，用于 transanno liftvcf 输入。

    contig 不在参考 FASTA 里的位点不写进 VCF（写了也只会被 transanno 判
    UNEXPECTED_REF，与真实比对失败混在一起），单独输出到 missing_out
    （默认 <id_file>.missing_contigs.tsv）；全部位点缺失时直接报错退出。
    """
    id_vcf_file = Path(f"{id_file}.vcf")
    if missing_out is None:
        missing_out = Path(f"{id_file}.missing_contigs.tsv")

    ref_fasta = Fasta(str(ref_fa_path))
    missing_contig_sites: list[tuple[str, str, int]] = []
    total_sites = 0
    with open(id_file) as id_inf, open(id_vcf_file, "w") as vcf_out:
        vcf_out.write("##fileformat=VCFv4.2\n")
        vcf_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for line in id_inf:
            each_id = line.strip()
            if not each_id:
                continue
            chrom, pos = parse_id_to_chrom_pos(each_id)
            total_sites += 1
            if chrom not in ref_fasta:
                missing_contig_sites.append((each_id, chrom, pos))
                continue
            ref_seq = fetch_ref_nucleotide(ref_fasta, chrom, pos)
            vcf_out.write(f"{chrom}\t{pos}\t{each_id}\t{ref_seq}\t.\t.\t.\t.\n")

    if missing_contig_sites:
        with open(missing_out, "w") as missing_out_f:
            for each_id, chrom, pos in missing_contig_sites:
                missing_out_f.write(f"{each_id}\t{chrom}\t{pos}\n")
        logger.warning(
            f"{len(missing_contig_sites)}/{total_sites} 个位点的 contig 不在参考 "
            f"{ref_fa_path} 中，未写入 VCF，已单独列出: {missing_out}；"
            "请检查 ID 的 contig 命名是否与参考一致"
        )
        if len(missing_contig_sites) == total_sites:
            raise ValueError(
                f"所有 {total_sites} 个位点的 contig 都不在参考 {ref_fa_path} 中，"
                "疑似 contig 命名整体不匹配，终止运行"
            )
    else:
        missing_out.unlink(missing_ok=True)

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
    expanded_end = df["end"] + flank
    chrom_limits = df["chrom"].astype(str).map(chrom_sizes)
    df["end"] = pd.concat(
        [expanded_end, chrom_limits.fillna(expanded_end)],
        axis=1,
    ).min(axis=1).astype(int)

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


def write_bed(df: pd.DataFrame, out_bed: Path, columns: Sequence[str]) -> None:
    """按 BED 三列格式写入文件。"""
    df.to_csv(out_bed, sep="\t", index=False, header=False, columns=columns)


def load_split_bed(split_bed: Path) -> pd.DataFrame:
    """读取并校验 split.bed 文件。"""
    if not split_bed.exists():
        raise FileNotFoundError(f"找不到 split.bed: {split_bed}")

    try:
        split_bed_df = pd.read_table(
            split_bed,
            header=None,
            names=["chrom", "split_start", "split_end", "new_chrom"],
            usecols=[0, 1, 2, 3],
        )
    except Exception as exc:  # pragma: no cover - pandas exception details vary
        raise ValueError(
            "split.bed 格式错误，必须包含 4 列: chrom, split_start, split_end, new_chrom"
        ) from exc

    split_bed_df["chrom"] = split_bed_df["chrom"].astype(str)
    split_bed_df["new_chrom"] = split_bed_df["new_chrom"].astype(str)
    split_bed_df["split_start"] = pd.to_numeric(
        split_bed_df["split_start"], errors="raise"
    ).astype(int)
    split_bed_df["split_end"] = pd.to_numeric(
        split_bed_df["split_end"], errors="raise"
    ).astype(int)
    return split_bed_df


def load_split_chrom_order(split_genome_fai: Path) -> list[str]:
    """读取 split 基因组染色体顺序。"""
    if not split_genome_fai.exists():
        raise FileNotFoundError(f"找不到 split genome fai: {split_genome_fai}")

    chrom_df = pd.read_table(split_genome_fai, header=None, names=["chrom"], usecols=[0])
    chrom_order = chrom_df["chrom"].astype(str).tolist()
    if not chrom_order:
        raise ValueError(f"split genome fai 为空: {split_genome_fai}")
    return chrom_order


def split_bed_dataframe(
    bed_df: pd.DataFrame,
    split_bed_df: pd.DataFrame,
    start_col: str,
    end_col: str,
) -> pd.DataFrame:
    """将 query genome BED 坐标映射到 split genome 坐标，正确处理跨边界重叠区间并进行 clamp。"""
    source_df = bed_df.copy()
    source_df["chrom"] = source_df["chrom"].astype(str)
    source_df["start"] = pd.to_numeric(source_df[start_col], errors="raise").astype(int)
    source_df["end"] = pd.to_numeric(source_df[end_col], errors="raise").astype(int)

    merge_df = source_df[["chrom", "start", "end"]].merge(split_bed_df, on="chrom")

    merge_df["overlap_start"] = merge_df[["start", "split_start"]].max(axis=1)
    merge_df["overlap_end"] = merge_df[["end", "split_end"]].min(axis=1)

    merge_df = merge_df[merge_df["overlap_start"] < merge_df["overlap_end"]].copy()

    merge_df["new_start"] = (merge_df["overlap_start"] - merge_df["split_start"]).astype(int)
    merge_df["new_end"] = (merge_df["overlap_end"] - merge_df["split_start"]).astype(int)
    return merge_df[["new_chrom", "new_start", "new_end"]]


def sort_split_output_by_fai(
    split_df: pd.DataFrame,
    split_chrom_order: list[str],
) -> pd.DataFrame:
    """按 split genome FAI 的染色体顺序排序 split 输出。"""
    if split_df.empty:
        return split_df

    ordered_df = split_df.copy()
    ordered_df["new_chrom"] = ordered_df["new_chrom"].astype(str)
    missing = sorted(set(ordered_df["new_chrom"].unique()) - set(split_chrom_order))
    if missing:
        raise ValueError(
            f"split 输出中存在不在 split genome fai 的染色体: {', '.join(missing)}"
        )

    ordered_df["new_chrom"] = pd.Categorical(
        ordered_df["new_chrom"],
        categories=split_chrom_order,
        ordered=True,
    )
    ordered_df = ordered_df.sort_values(
        by=["new_chrom", "new_start", "new_end"],
        kind="mergesort",
    ).reset_index(drop=True)
    ordered_df["new_chrom"] = ordered_df["new_chrom"].astype(str)
    return ordered_df


def resolve_split_genome_fai(
    split_bed: Path,
    split_genome_fai: Path | None,
) -> Path:
    """解析 split genome fai；未显式提供时使用 split.bed 同目录的 genome.fa.fai。"""
    if split_genome_fai is not None:
        return split_genome_fai

    inferred_fai = split_bed.parent / "genome.fa.fai"
    if inferred_fai.exists():
        return inferred_fai

    raise ValueError(
        "启用 --split-bed 时需要 --split-genome-fai，"
        "或在 split.bed 同目录提供 genome.fa.fai"
    )


def write_liftover_outputs(
    sorted_bed: pd.DataFrame,
    snpcalling_sorted: pd.DataFrame,
    outdir: Path,
    probe_name: str,
    split_bed: Path | None,
    split_genome_fai: Path | None,
) -> None:
    """写入 liftover ID、坐标表和 BED；启用 split 时仅 BED 坐标转换到 split 基因组。"""
    probe_id_file = outdir / f"{probe_name}.id"
    sorted_bed.to_csv(probe_id_file, sep="\t", index=False, header=False, columns=["pos_id"])
    logger.info(f"已写入 ID 文件: {probe_id_file}")

    probe_pos_file = outdir / f"{probe_name}.pos.tsv"
    sorted_bed.to_csv(
        probe_pos_file,
        sep="\t",
        index=False,
        header=False,
        columns=["chrom", "pos", "id"],
    )
    logger.info(f"已写入位点坐标文件: {probe_pos_file}")

    probe_bed = outdir / f"{probe_name}.bed"
    snpcalling_bed = outdir / f"{probe_name}.snpcalling.bed"
    if split_bed is None:
        write_bed(sorted_bed, probe_bed, columns=["chrom", "start", "pos"])
        logger.info(f"已写入 BED 文件: {probe_bed}")
        write_bed(snpcalling_sorted, snpcalling_bed, columns=["chrom", "start", "end"])
        logger.info(f"已写入 snpcalling BED: {snpcalling_bed}")
        return

    split_bed_df = load_split_bed(split_bed)
    resolved_split_genome_fai = resolve_split_genome_fai(split_bed, split_genome_fai)
    split_chrom_order = load_split_chrom_order(resolved_split_genome_fai)

    split_target_df = split_bed_dataframe(
        sorted_bed,
        split_bed_df,
        start_col="start",
        end_col="pos",
    )
    split_target_df = sort_split_output_by_fai(split_target_df, split_chrom_order)
    split_snpcalling_df = split_bed_dataframe(
        snpcalling_sorted,
        split_bed_df,
        start_col="start",
        end_col="end",
    )
    split_snpcalling_df = sort_split_output_by_fai(
        split_snpcalling_df,
        split_chrom_order,
    )

    write_bed(split_target_df, probe_bed, columns=["new_chrom", "new_start", "new_end"])
    logger.info(f"已写入 split BED 文件: {probe_bed}")
    write_bed(
        split_snpcalling_df,
        snpcalling_bed,
        columns=["new_chrom", "new_start", "new_end"],
    )
    logger.info(f"已写入 split snpcalling BED: {snpcalling_bed}")


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
    split_bed: Annotated[
        Path | None,
        typer.Option(help="split.bed 文件路径；提供后 BED 输出 split 坐标"),
    ] = None,
    split_genome_fai: Annotated[
        Path | None,
        typer.Option(
            help="split 基因组 fai（.fai），用于按染色体顺序排序 split 输出；"
            "未提供时会尝试使用 split.bed 同目录下的 genome.fa.fai"
        ),
    ] = None,
    probe_name: Annotated[
        str | None,
        typer.Option(
            help="输出文件名前缀；默认使用 id_file 的 stem。"
            "当 id_file 被 stage 为固定名时由流程传入原始 stem。"
        ),
    ] = None,
) -> None:
    """根据 ID 文件进行 liftover 并生成 BED/ID 文件。"""
    if flank < 0:
        raise typer.BadParameter(f"flank 必须大于等于 0，当前值: {flank}")

    outdir.mkdir(exist_ok=True, parents=True)

    # 1. 生成 VCF；contig 缺失的位点单独列到 outdir，随 out/* 一起发布
    vcf = make_id_vcf(
        id_file,
        ref_fa,
        missing_out=outdir / f"{id_file.name}.missing_contigs.tsv",
    )

    # 2. 运行 transanno liftvcf
    probe = probe_name if probe_name else id_file.stem
    lift_over_vcf = outdir / f"liftover.{probe}.id.vcf.gz"
    rejected_vcf = outdir / f"rejected.{probe}.id.vcf.gz"

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

    # 5. 生成 snpcalling BED (slop + merge)
    chrom_sizes_dict = dict(zip(chrom_df["chrom"], chrom_df["chrom_size"], strict=True))
    snpcalling_df = slop_and_merge(sorted_bed, chrom_sizes_dict, flank)
    snpcalling_sorted = sort_bed_by_fai(snpcalling_df, chrom_df)

    # 6. 写入输出；启用 split 时仅 BED 坐标转换到 split 基因组。
    write_liftover_outputs(
        sorted_bed=sorted_bed,
        snpcalling_sorted=snpcalling_sorted,
        outdir=outdir,
        probe_name=probe,
        split_bed=split_bed,
        split_genome_fai=split_genome_fai,
    )


if __name__ == "__main__":
    app()
