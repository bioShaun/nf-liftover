#!/usr/bin/env python3
"""Merge transanno liftover results with strict BLAST-rescued failures."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd
from loguru import logger
from pyfaidx import Fasta

from liftover_by_id import (
    fetch_ref_nucleotide,
    load_chrom_sizes_from_fai,
    slop_and_merge,
    sort_bed_by_fai,
    write_liftover_outputs,
)

REASON_RE = re.compile(r"(?:^|[;.])FAILED_REASON=([^;]+)")


def prefixed_path(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}{suffix}")


def panel_path(outdir: Path, panel: str, suffix: str) -> Path:
    return outdir / f"{panel}{suffix}"


def read_input_ids(id_file: Path) -> list[str]:
    return [line.strip() for line in id_file.read_text(encoding="utf-8").splitlines() if line.strip()]


def read_transanno_pos(pos_tsv: Path) -> pd.DataFrame:
    df = pd.read_table(pos_tsv, header=None, names=["chrom", "pos", "id"], dtype={"chrom": str, "id": str})
    df["pos"] = pd.to_numeric(df["pos"], errors="raise").astype(int)
    df["method"] = "transanno"
    return df[["chrom", "pos", "id", "method"]]


def read_mapping_targets(mapping_tsv: Path) -> set[str]:
    df = pd.read_table(mapping_tsv, header=None, usecols=[1], names=["target"], dtype=str)
    return set(df["target"].dropna().astype(str))


def read_rejected_reasons(rejected_vcf: Path | None) -> dict[str, str]:
    if rejected_vcf is None or not rejected_vcf.exists():
        return {}
    opener = gzip.open if rejected_vcf.suffix == ".gz" else open
    reasons: dict[str, str] = {}
    with opener(rejected_vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            marker_id = fields[2]
            match = REASON_RE.search(fields[7])
            reasons[marker_id] = match.group(1) if match else "TRANSANNO_REJECTED"
    return reasons


def read_rejected_alleles(rejected_vcf: Path | None) -> dict[str, tuple[str, str | None]]:
    """从 rejected VCF 取源位点 REF/ALT（ALT 为 '.' 时视为无 ALT）。"""
    if rejected_vcf is None or not rejected_vcf.exists():
        return {}
    opener = gzip.open if rejected_vcf.suffix == ".gz" else open
    alleles: dict[str, tuple[str, str | None]] = {}
    with opener(rejected_vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            ref = fields[3].upper()
            alt = fields[4].split(",")[0].upper()
            alleles[fields[2]] = (ref, alt if alt not in {"", ".", "-"} else None)
    return alleles


SELECTION_NUMERIC_COLUMNS = [
    "rank",
    "pos",
    "bitscore",
    "match_ratio",
    "align_len",
    "mismatches",
    "pident",
    "global_best_bitscore",
    "global_best_pident",
    "global_second_best_bitscore",
]


def read_selection(selection_tsv: Path) -> pd.DataFrame:
    if not selection_tsv.exists() or selection_tsv.stat().st_size == 0:
        return pd.DataFrame()
    df = pd.read_table(selection_tsv, dtype={"id": str, "chrom": str})
    for column in SELECTION_NUMERIC_COLUMNS:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="raise")
    if "pident" not in df.columns and {"mismatches", "align_len"} <= set(df.columns):
        # 过渡期近似：旧格式 selection.tsv 没有 pident 列
        align_len = df["align_len"].where(df["align_len"] > 0)
        df["pident"] = 100.0 * (1.0 - df["mismatches"] / align_len)
        logger.warning("selection.tsv 缺少 pident 列，用 1 - mismatches/align_len 近似（过渡期行为）")
    return df


def _ambiguity_scores(group: pd.DataFrame, top_row: pd.Series) -> tuple[float, float] | None:
    """返回同一候选空间内的最佳和次优 bitscore。"""
    global_best = top_row.get("global_best_bitscore")
    global_second = top_row.get("global_second_best_bitscore")
    if pd.notna(global_best) and pd.notna(global_second):
        return float(global_best), float(global_second)
    others = group.drop(index=top_row.name)
    if others.empty:
        return None
    return float(top_row["bitscore"]), float(others.sort_values("rank").iloc[0]["bitscore"])


def _is_ambiguous_top(
    group: pd.DataFrame,
    top_row: pd.Series,
    marker_id: str,
    *,
    tied_margin: float,
    tied_margin_frac: float,
) -> bool:
    """margin 阈值判歧义：best - second < max(tied_margin, tied_margin_frac * best)。

    margin ≈ 0 的位点（多倍体 A/B 拷贝）是真歧义，判 True 交由下游标注处理。
    """
    scores = _ambiguity_scores(group, top_row)
    if scores is None:
        logger.warning(
            f"{marker_id}: group 内只有 1 条记录且无完整 global best/second bitscore，"
            "无法判断 top-hit 歧义，按无歧义放行"
        )
        return False
    best, second = scores
    return best - second < max(tied_margin, tied_margin_frac * best)


def _parse_alleles_field(value: object) -> tuple[str, str | None] | None:
    """解析 selection.tsv 的 alleles 列（如 ``A/G``）；``-/-`` 等无效值返回 None。"""
    text = str(value or "").strip().upper()
    if "/" not in text:
        return None
    ref, alt = text.split("/", 1)
    if ref in {"", "-", "."}:
        return None
    return ref, alt if alt not in {"", "-", "."} else None


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def classify_blast_rescue(
    selection_df: pd.DataFrame,
    failed_ids: set[str],
    mapping_targets: set[str],
    *,
    min_query_coverage: float = 0.95,
    min_pident: float = 95.0,
    max_constraint_pident_drop: float = 2.0,
    tied_margin: float = 20.0,
    tied_margin_frac: float = 0.05,
    target_fasta: Fasta | None = None,
    source_alleles: dict[str, tuple[str, str | None]] | None = None,
) -> tuple[pd.DataFrame, dict[str, str], dict[str, int]]:
    accepted = []
    blast_reasons = {marker_id: "NO_BLAST_RESCUE_HIT" for marker_id in failed_ids}
    ref_stats = {"ref_check_passed": 0, "ref_check_ref_is_alt": 0, "ref_check_mismatch": 0, "ref_check_skipped": 0}
    if selection_df.empty:
        return pd.DataFrame(columns=["chrom", "pos", "id", "method"]), blast_reasons, ref_stats

    has_pident = "pident" in selection_df.columns
    has_global_best = "global_best_pident" in selection_df.columns
    if not has_pident:
        logger.warning("selection.tsv 无 pident 且无法近似，BLAST_LOW_IDENTITY 闸门降级跳过（旧格式）")
    if not has_global_best:
        logger.warning("selection.tsv 缺少 global_best_* 列，BLAST_CONSTRAINT_OVERRODE_BEST_HIT 闸门降级跳过（旧格式）")

    for marker_id, group in selection_df[selection_df["id"].isin(failed_ids)].groupby("id", sort=False):
        top = group[group["rank"] == 1]
        if top.empty:
            blast_reasons[marker_id] = "NO_RANK1_BLAST_HIT"
            continue
        row = top.iloc[0]
        checks = []
        if float(row["match_ratio"]) < min_query_coverage:
            checks.append("BLAST_LOW_MATCH_RATIO")
        if has_pident and float(row["pident"]) < min_pident:
            checks.append("BLAST_LOW_IDENTITY")
        if has_global_best and has_pident:
            global_best_pident = row.get("global_best_pident")
            if pd.notna(global_best_pident) and float(global_best_pident) - float(row["pident"]) > max_constraint_pident_drop:
                checks.append("BLAST_CONSTRAINT_OVERRODE_BEST_HIT")
        # matched_but_not_global_best 视为"匹配"，是否采信交给上面的 pident 落差闸门判断
        if str(row.get("chr_map_status", "")) not in ("matched", "matched_but_not_global_best"):
            checks.append("BLAST_CHR_MAP_NOT_MATCHED")
        if str(row["chrom"]) not in mapping_targets:
            checks.append("BLAST_TARGET_CHROM_NOT_IN_MAPPING")
        if _is_ambiguous_top(group, row, marker_id, tied_margin=tied_margin, tied_margin_frac=tied_margin_frac):
            checks.append("BLAST_AMBIGUOUS_TOP_HIT")

        if checks:
            blast_reasons[marker_id] = "+".join(checks)
            continue

        # REF 核对：只对通过上述闸门的待采信 rescue 做（与 transanno 对齐的验收强度）
        if target_fasta is not None:
            alleles = _parse_alleles_field(row.get("alleles")) or (source_alleles or {}).get(marker_id)
            if alleles is None:
                ref_stats["ref_check_skipped"] += 1
                logger.warning(f"{marker_id}: 无 REF/ALT 信息（selection alleles 与 rejected VCF 均缺），跳过 REF 核对")
            else:
                ref, alt = alleles
                if str(row.get("strand", "+")) == "-":
                    ref = _reverse_complement(ref)
                    alt = _reverse_complement(alt) if alt is not None else None
                target_base = fetch_ref_nucleotide(target_fasta, str(row["chrom"]), int(row["pos"]))
                if target_base == ref:
                    ref_stats["ref_check_passed"] += 1
                elif alt is not None and target_base == alt:
                    ref_stats["ref_check_ref_is_alt"] += 1
                    blast_reasons[marker_id] = "BLAST_REF_IS_ALT"
                    continue
                else:
                    ref_stats["ref_check_mismatch"] += 1
                    blast_reasons[marker_id] = "BLAST_REF_MISMATCH"
                    continue

        accepted.append({"chrom": str(row["chrom"]), "pos": int(row["pos"]), "id": marker_id, "method": "blast_rescue"})
        blast_reasons[marker_id] = "BLAST_RESCUED"

    return pd.DataFrame(accepted, columns=["chrom", "pos", "id", "method"]), blast_reasons, ref_stats


def make_status_table(
    input_ids: list[str],
    transanno_df: pd.DataFrame,
    rescue_df: pd.DataFrame,
    transanno_reasons: dict[str, str],
    blast_reasons: dict[str, str],
) -> pd.DataFrame:
    transanno_ids = set(transanno_df["id"].astype(str))
    rescue_ids = set(rescue_df["id"].astype(str))
    rows = []
    for marker_id in input_ids:
        if marker_id in transanno_ids:
            rows.append({"id": marker_id, "status": "transanno_success", "failed_reason": ""})
        elif marker_id in rescue_ids:
            rows.append({"id": marker_id, "status": "blast_rescue", "failed_reason": transanno_reasons.get(marker_id, "TRANSANNO_FAILED")})
        else:
            t_reason = transanno_reasons.get(marker_id, "TRANSANNO_FAILED")
            b_reason = blast_reasons.get(marker_id, "NO_BLAST_RESCUE_HIT")
            rows.append({"id": marker_id, "status": "failed", "failed_reason": f"{t_reason};{b_reason}"})
    return pd.DataFrame(rows)


def write_summary(status_df: pd.DataFrame, outdir: Path, panel: str, ref_stats: dict[str, int] | None = None) -> None:
    total = len(status_df)
    counts = status_df["status"].value_counts().to_dict()
    rows = []
    for status in ["transanno_success", "blast_rescue", "failed"]:
        count = int(counts.get(status, 0))
        rows.append({"metric": status, "count": count, "percent": count / total if total else 0})
    if ref_stats:
        checked = ref_stats["ref_check_passed"] + ref_stats["ref_check_ref_is_alt"] + ref_stats["ref_check_mismatch"]
        for metric in ["ref_check_passed", "ref_check_ref_is_alt", "ref_check_mismatch", "ref_check_skipped"]:
            count = ref_stats[metric]
            rows.append({"metric": metric, "count": count, "percent": count / total if total else 0})
        rows.append({
            "metric": "ref_check_pass_rate",
            "count": ref_stats["ref_check_passed"],
            "percent": ref_stats["ref_check_passed"] / checked if checked else 0,
        })
    pd.DataFrame(rows).to_csv(panel_path(outdir, panel, ".summary.tsv"), sep="\t", index=False)

    failed = status_df[status_df["status"] == "failed"]
    if failed.empty:
        reason_df = pd.DataFrame(columns=["failed_reason", "count"])
    else:
        reason_df = failed["failed_reason"].value_counts().rename_axis("failed_reason").reset_index(name="count")
    reason_df.to_csv(panel_path(outdir, panel, ".failed-reasons.tsv"), sep="\t", index=False)
    status_df.to_csv(panel_path(outdir, panel, ".status.tsv"), sep="\t", index=False)


def write_combined_outputs(
    combined_df: pd.DataFrame,
    query_fai: Path,
    outdir: Path,
    panel: str,
    *,
    flank: int,
    split_bed: Path | None = None,
    split_genome_fai: Path | None = None,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    chrom_df = load_chrom_sizes_from_fai(query_fai)
    bed_df = combined_df.copy()
    bed_df["start"] = bed_df["pos"] - 1
    bed_df["pos_id"] = bed_df["chrom"].astype(str) + "_" + bed_df["pos"].astype(str)
    sorted_df = sort_bed_by_fai(bed_df, chrom_df)

    chrom_sizes = dict(zip(chrom_df["chrom"], chrom_df["chrom_size"], strict=True))
    snpcalling = slop_and_merge(sorted_df, chrom_sizes, flank)
    snpcalling_sorted = sort_bed_by_fai(snpcalling, chrom_df)

    write_liftover_outputs(
        sorted_bed=sorted_df,
        snpcalling_sorted=snpcalling_sorted,
        outdir=outdir,
        probe_name=panel,
        split_bed=split_bed,
        split_genome_fai=split_genome_fai,
    )
    sorted_df.to_csv(panel_path(outdir, panel, ".method.tsv"), sep="\t", index=False, columns=["chrom", "pos", "id", "method"])


def merge_blast_rescue(
    *,
    id_file: Path,
    transanno_pos: Path,
    rejected_vcf: Path | None,
    blast_selection: Path,
    mapping_tsv: Path,
    query_fai: Path,
    outdir: Path,
    panel: str,
    flank: int = 100,
    min_query_coverage: float = 0.95,
    min_pident: float = 95.0,
    max_constraint_pident_drop: float = 2.0,
    tied_margin: float = 20.0,
    tied_margin_frac: float = 0.05,
    target_fasta: Path | None = None,
    split_bed: Path | None = None,
    split_genome_fai: Path | None = None,
) -> None:
    input_ids = read_input_ids(id_file)
    transanno_df = read_transanno_pos(transanno_pos)
    transanno_df = transanno_df[transanno_df["id"].isin(input_ids)].copy()
    transanno_ids = set(transanno_df["id"].astype(str))
    failed_ids = set(input_ids) - transanno_ids

    rescue_df, blast_reasons, ref_stats = classify_blast_rescue(
        read_selection(blast_selection),
        failed_ids,
        read_mapping_targets(mapping_tsv),
        min_query_coverage=min_query_coverage,
        min_pident=min_pident,
        max_constraint_pident_drop=max_constraint_pident_drop,
        tied_margin=tied_margin,
        tied_margin_frac=tied_margin_frac,
        target_fasta=Fasta(str(target_fasta)) if target_fasta is not None else None,
        source_alleles=read_rejected_alleles(rejected_vcf),
    )
    transanno_reasons = read_rejected_reasons(rejected_vcf)
    status_df = make_status_table(input_ids, transanno_df, rescue_df, transanno_reasons, blast_reasons)

    combined_df = pd.concat([transanno_df, rescue_df], ignore_index=True)
    write_combined_outputs(
        combined_df,
        query_fai,
        outdir,
        panel,
        flank=flank,
        split_bed=split_bed,
        split_genome_fai=split_genome_fai,
    )
    write_summary(status_df, outdir, panel, ref_stats=ref_stats if target_fasta is not None else None)
    rescue_df.to_csv(panel_path(outdir, panel, ".blast-rescue.accepted.tsv"), sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--id", required=True, type=Path, help="Original input .id file")
    parser.add_argument("--transanno-pos", required=True, type=Path, help="Transanno liftover <panel>.pos.tsv")
    parser.add_argument("--rejected-vcf", type=Path, help="Transanno rejected.<panel>.id.vcf.gz")
    parser.add_argument("--blast-selection", required=True, type=Path, help="realign_blast *.selection.tsv")
    parser.add_argument("--mapping", required=True, type=Path, help="source-to-target mapping TSV")
    parser.add_argument("--query-fai", required=True, type=Path, help="Target FASTA .fai")
    parser.add_argument("--outdir", type=Path, help="Output directory; filenames stay <panel>.id/.bed/.snpcalling.bed")
    parser.add_argument("--panel", help="Output filename stem; defaults to input .id stem")
    parser.add_argument("--out-prefix", type=Path, help="Deprecated: use --outdir and --panel")
    parser.add_argument("--flank", type=int, default=100, help="snpcalling BED flank size")
    parser.add_argument(
        "--min-query-coverage",
        "--min-match-ratio",
        dest="min_query_coverage",
        type=float,
        default=0.95,
        help="Minimum BLAST rescue query coverage (align_len / informative_len); --min-match-ratio kept as deprecated alias",
    )
    parser.add_argument("--min-pident", type=float, default=95.0, help="Minimum BLAST rescue percent identity")
    parser.add_argument(
        "--max-constraint-pident-drop",
        type=float,
        default=2.0,
        help="Reject when the constrained hit's pident trails global_best_pident by more than this (percentage points)",
    )
    parser.add_argument("--tied-margin", type=float, default=20.0, help="Top-hit ambiguity margin: absolute bitscore floor")
    parser.add_argument("--tied-margin-frac", type=float, default=0.05, help="Top-hit ambiguity margin: fraction of best bitscore")
    parser.add_argument("--target-fasta", required=True, type=Path, help="Target genome FASTA; rescued sites are REF-verified against it")
    parser.add_argument("--split-bed", type=Path, help="Optional split.bed for split-coordinate BED output")
    parser.add_argument("--split-genome-fai", type=Path, help="Optional split genome FASTA index for split BED sorting")
    args = parser.parse_args()

    if args.outdir is None:
        if args.out_prefix is None:
            parser.error("--outdir is required")
        # Backwards-compatible shim; new workflows should use a separate outdir.
        outdir = args.out_prefix.parent
        panel = args.out_prefix.name
    else:
        outdir = args.outdir
        panel = args.panel or args.id.stem

    merge_blast_rescue(
        id_file=args.id,
        transanno_pos=args.transanno_pos,
        rejected_vcf=args.rejected_vcf,
        blast_selection=args.blast_selection,
        mapping_tsv=args.mapping,
        query_fai=args.query_fai,
        outdir=outdir,
        panel=panel,
        flank=args.flank,
        min_query_coverage=args.min_query_coverage,
        min_pident=args.min_pident,
        max_constraint_pident_drop=args.max_constraint_pident_drop,
        tied_margin=args.tied_margin,
        tied_margin_frac=args.tied_margin_frac,
        target_fasta=args.target_fasta,
        split_bed=args.split_bed,
        split_genome_fai=args.split_genome_fai,
    )


if __name__ == "__main__":
    main()
