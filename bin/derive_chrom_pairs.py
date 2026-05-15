#!/usr/bin/env python3
"""Derive query-to-target chromosome pairs from FASTA index files."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class ChromRecord:
    name: str
    length: int


def read_fai(path: Path) -> list[ChromRecord]:
    records: list[ChromRecord] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_no}: expected at least 2 tab-separated columns")
            records.append(ChromRecord(fields[0], int(fields[1])))
    if not records:
        raise ValueError(f"No chromosomes found in {path}")
    return records


def read_mapping(path: Path) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(f"{path}:{line_no}: mapping must contain exactly two columns")
            pairs.append((fields[0], fields[1]))
    if not pairs:
        raise ValueError(f"No chromosome pairs found in {path}")
    return pairs


def numeric_suffix(name: str) -> str | None:
    match = re.search(r"(\d+)$", name)
    if not match:
        return None
    return str(int(match.group(1)))


def derive_by_order(ref_records: list[ChromRecord], query_records: list[ChromRecord]) -> list[tuple[str, str]]:
    if len(ref_records) != len(query_records):
        raise ValueError(
            f"Cannot derive pairs by order: ref has {len(ref_records)} chromosomes, "
            f"query has {len(query_records)} chromosomes"
        )
    return [(ref.name, query.name) for ref, query in zip(ref_records, query_records, strict=True)]


def derive_by_suffix(ref_records: list[ChromRecord], query_records: list[ChromRecord]) -> list[tuple[str, str]]:
    query_by_suffix: dict[str, ChromRecord] = {}
    for record in query_records:
        suffix = numeric_suffix(record.name)
        if suffix is None:
            continue
        if suffix in query_by_suffix:
            raise ValueError(f"Ambiguous query chromosome suffix {suffix!r}")
        query_by_suffix[suffix] = record

    pairs: list[tuple[str, str]] = []
    missing: list[str] = []
    for ref in ref_records:
        suffix = numeric_suffix(ref.name)
        if suffix is None or suffix not in query_by_suffix:
            missing.append(ref.name)
            continue
        pairs.append((ref.name, query_by_suffix[suffix].name))

    if missing:
        raise ValueError(f"Cannot match query chromosomes by numeric suffix for: {', '.join(missing)}")
    return pairs


def validate_pairs(
    pairs: Iterable[tuple[str, str]],
    ref_records: list[ChromRecord],
    query_records: list[ChromRecord],
) -> list[tuple[str, str]]:
    ref_names = {record.name for record in ref_records}
    query_names = {record.name for record in query_records}
    validated: list[tuple[str, str]] = []

    for ref_name, query_name in pairs:
        if ref_name not in ref_names:
            raise ValueError(f"Mapping references unknown ref chromosome: {ref_name}")
        if query_name not in query_names:
            raise ValueError(f"Mapping references unknown query chromosome: {query_name}")
        validated.append((ref_name, query_name))

    return validated


def write_pairs(path: Path, pairs: Iterable[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for ref_name, query_name in pairs:
            handle.write(f"{ref_name}\t{query_name}\n")


def derive_chrom_pairs(
    ref_fai: Path,
    query_fai: Path,
    output: Path,
    *,
    mapping: Path | None = None,
    strategy: str = "order",
) -> None:
    ref_records = read_fai(ref_fai)
    query_records = read_fai(query_fai)

    if mapping is not None:
        pairs = read_mapping(mapping)
    elif strategy == "order":
        pairs = derive_by_order(ref_records, query_records)
    elif strategy == "suffix":
        pairs = derive_by_suffix(ref_records, query_records)
    else:
        raise ValueError(f"Unknown pairing strategy: {strategy}")

    write_pairs(output, validate_pairs(pairs, ref_records, query_records))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-fai", required=True, type=Path, help="Reference/query-origin FASTA index")
    parser.add_argument("--query-fai", required=True, type=Path, help="Target/new FASTA index")
    parser.add_argument("--output", required=True, type=Path, help="Output chromosome pairs TSV")
    parser.add_argument("--mapping", type=Path, help="Optional user-provided two-column TSV mapping")
    parser.add_argument(
        "--strategy",
        choices=("order", "suffix"),
        default="order",
        help="Automatic derivation strategy when --mapping is absent",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    derive_chrom_pairs(
        args.ref_fai,
        args.query_fai,
        args.output,
        mapping=args.mapping,
        strategy=args.strategy,
    )


if __name__ == "__main__":
    main()
