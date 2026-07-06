#!/usr/bin/env python3
"""Merge transanno liftover results with strict BLAST-rescued failures."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd

from liftover_by_id import (
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


def read_selection(selection_tsv: Path) -> pd.DataFrame:
    if not selection_tsv.exists() or selection_tsv.stat().st_size == 0:
        return pd.DataFrame()
    df = pd.read_table(selection_tsv, dtype={"id": str, "chrom": str})
    for column in ["rank", "pos", "bitscore", "match_ratio"]:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="raise")
    return df


def _has_tied_top(group: pd.DataFrame) -> bool:
    top = group[group["rank"] == 1]
    if top.empty:
        return False
    top_bitscore = float(top.iloc[0]["bitscore"])
    return (group["bitscore"].astype(float) == top_bitscore).sum() > 1


def classify_blast_rescue(
    selection_df: pd.DataFrame,
    failed_ids: set[str],
    mapping_targets: set[str],
    *,
    min_match_ratio: float = 0.95,
) -> tuple[pd.DataFrame, dict[str, str]]:
    accepted = []
    blast_reasons = {marker_id: "NO_BLAST_RESCUE_HIT" for marker_id in failed_ids}
    if selection_df.empty:
        return pd.DataFrame(columns=["chrom", "pos", "id", "method"]), blast_reasons

    for marker_id, group in selection_df[selection_df["id"].isin(failed_ids)].groupby("id", sort=False):
        top = group[group["rank"] == 1]
        if top.empty:
            blast_reasons[marker_id] = "NO_RANK1_BLAST_HIT"
            continue
        row = top.iloc[0]
        checks = []
        if float(row["match_ratio"]) < min_match_ratio:
            checks.append("BLAST_LOW_MATCH_RATIO")
        if str(row.get("chr_map_status", "")) != "matched":
            checks.append("BLAST_CHR_MAP_NOT_MATCHED")
        if str(row["chrom"]) not in mapping_targets:
            checks.append("BLAST_TARGET_CHROM_NOT_IN_MAPPING")
        if _has_tied_top(group):
            checks.append("BLAST_TIED_TOP_HIT")

        if checks:
            blast_reasons[marker_id] = "+".join(checks)
            continue

        accepted.append({"chrom": str(row["chrom"]), "pos": int(row["pos"]), "id": marker_id, "method": "blast_rescue"})
        blast_reasons[marker_id] = "BLAST_RESCUED"

    return pd.DataFrame(accepted, columns=["chrom", "pos", "id", "method"]), blast_reasons


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


def write_summary(status_df: pd.DataFrame, outdir: Path, panel: str) -> None:
    total = len(status_df)
    counts = status_df["status"].value_counts().to_dict()
    rows = []
    for status in ["transanno_success", "blast_rescue", "failed"]:
        count = int(counts.get(status, 0))
        rows.append({"metric": status, "count": count, "percent": count / total if total else 0})
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
    min_match_ratio: float = 0.95,
    split_bed: Path | None = None,
    split_genome_fai: Path | None = None,
) -> None:
    input_ids = read_input_ids(id_file)
    transanno_df = read_transanno_pos(transanno_pos)
    transanno_df = transanno_df[transanno_df["id"].isin(input_ids)].copy()
    transanno_ids = set(transanno_df["id"].astype(str))
    failed_ids = set(input_ids) - transanno_ids

    rescue_df, blast_reasons = classify_blast_rescue(
        read_selection(blast_selection),
        failed_ids,
        read_mapping_targets(mapping_tsv),
        min_match_ratio=min_match_ratio,
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
    write_summary(status_df, outdir, panel)
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
    parser.add_argument("--min-match-ratio", type=float, default=0.95, help="Minimum BLAST rescue match_ratio")
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
        min_match_ratio=args.min_match_ratio,
        split_bed=args.split_bed,
        split_genome_fai=args.split_genome_fai,
    )


if __name__ == "__main__":
    main()
