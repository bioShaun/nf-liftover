#!/usr/bin/env python3
"""Restore original query coordinates for PAF records generated from seqkit sliding windows."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


SLIDING_RE = re.compile(r"^(?P<chrom>.+)_sliding:(?P<start>\d+)-(?P<end>\d+)$")


def load_fai_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_no}: expected at least two columns")
            lengths[fields[0]] = int(fields[1])
    if not lengths:
        raise ValueError(f"No sequence lengths found in {path}")
    return lengths


def parse_sliding_query_name(name: str) -> tuple[str, int]:
    match = SLIDING_RE.match(name)
    if not match:
        raise ValueError(f"Cannot parse sliding query name: {name}")
    chrom = match.group("chrom")
    offset = int(match.group("start")) - 1
    return chrom, offset


def restore_split_paf(input_paf: Path, query_fai: Path, output_paf: Path) -> None:
    lengths = load_fai_lengths(query_fai)
    output_paf.parent.mkdir(parents=True, exist_ok=True)

    with input_paf.open(encoding="utf-8") as inf, output_paf.open("w", encoding="utf-8") as out:
        for line_no, line in enumerate(inf, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError(f"{input_paf}:{line_no}: expected at least 12 PAF columns")

            chrom, offset = parse_sliding_query_name(fields[0])
            if chrom not in lengths:
                raise ValueError(f"{input_paf}:{line_no}: {chrom!r} not found in {query_fai}")

            fields[0] = chrom
            fields[1] = str(lengths[chrom])
            fields[2] = str(int(fields[2]) + offset)
            fields[3] = str(int(fields[3]) + offset)
            out.write("\t".join(fields) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_paf", type=Path, help="Input split-window PAF")
    parser.add_argument("query_fai", type=Path, help="Original query FASTA index")
    parser.add_argument("output_paf", type=Path, help="Output restored PAF")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    restore_split_paf(args.input_paf, args.query_fai, args.output_paf)


if __name__ == "__main__":
    main()
