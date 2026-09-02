#!/usr/bin/env python3
"""Repair missing VCF INFO and contig header declarations without changing records."""
from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import TextIO

INFO_ID_RE = re.compile(r"^##INFO=<ID=([^,>]+)")
CONTIG_ID_RE = re.compile(r"^##contig=<ID=([^,>]+)")


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="strict")
    return path.open("r", encoding="utf-8")


def repair_header(vcf: Path, fai: Path, output_header: Path) -> dict[str, object]:
    header: list[str] = []
    chrom_header: str | None = None
    existing_info: set[str] = set()
    existing_contigs: set[str] = set()
    observed: dict[str, set[str]] = defaultdict(set)
    records = 0

    with open_text(vcf) as handle:
        for line in handle:
            if line.startswith("##"):
                header.append(line.rstrip("\n"))
                match = INFO_ID_RE.match(line)
                if match:
                    existing_info.add(match.group(1))
                match = CONTIG_ID_RE.match(line)
                if match:
                    existing_contigs.add(match.group(1))
                continue
            if line.startswith("#CHROM"):
                chrom_header = line.rstrip("\n")
                continue
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record with {len(fields)} columns")
            records += 1
            info = fields[7]
            if info == ".":
                continue
            for item in info.split(";"):
                if not item:
                    continue
                if "=" in item:
                    key, _ = item.split("=", 1)
                    observed[key].add("value")
                else:
                    observed[item].add("flag")

    if chrom_header is None:
        raise ValueError("VCF #CHROM header is missing")
    mixed = {key: sorted(kinds) for key, kinds in observed.items() if len(kinds) > 1 and key not in existing_info}
    if mixed:
        raise ValueError(f"Undeclared INFO keys used as both flag and value: {mixed}")

    contig_lines: list[str] = []
    with fai.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Malformed FAI line: {line.rstrip()}")
            contig, length = fields[0], fields[1]
            if contig not in existing_contigs:
                contig_lines.append(f"##contig=<ID={contig},length={length}>")

    info_lines: list[str] = []
    for key in sorted(observed):
        if key in existing_info:
            continue
        if observed[key] == {"flag"}:
            info_lines.append(
                f'##INFO=<ID={key},Number=0,Type=Flag,Description="Header inferred from bare INFO usage; records unchanged">'
            )
        else:
            info_lines.append(
                f'##INFO=<ID={key},Number=.,Type=String,Description="Header inferred from INFO=value usage; records unchanged">'
            )

    output_header.parent.mkdir(parents=True, exist_ok=True)
    output_header.write_text(
        "\n".join(header + contig_lines + info_lines + [chrom_header]) + "\n",
        encoding="utf-8",
    )
    return {
        "records_scanned": records,
        "observed_info_keys": sorted(observed),
        "existing_info_keys": sorted(existing_info),
        "missing_info_added": len(info_lines),
        "missing_contigs_added": len(contig_lines),
        "output_header": str(output_header),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--fai", type=Path, required=True)
    parser.add_argument("--output-header", type=Path, required=True)
    parser.add_argument("--summary-json", type=Path)
    args = parser.parse_args()
    summary = repair_header(args.vcf, args.fai, args.output_header)
    text = json.dumps(summary, indent=2, sort_keys=True)
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
