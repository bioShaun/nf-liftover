#!/usr/bin/env python3
"""Chain metadata + validation helpers for nf-liftover."""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import sys
from pathlib import Path
from typing import Any

_SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


def sha256_file(path: Path | str) -> tuple[str, int]:
    h = hashlib.sha256()
    size = 0
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
            size += len(chunk)
    return h.hexdigest(), size


def load_fai(path: Path | str) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            raise ValueError(f"Invalid FAI line in {path}: {line!r}")
        lengths[parts[0]] = int(parts[1])
    return lengths


def parse_chain_headers(path: Path | str) -> list[dict[str, Any]]:
    """Parse UCSC chain headers: tName/tSize are source, qName/qSize are target."""
    blocks: list[dict[str, Any]] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.startswith("chain"):
                continue
            fields = line.split()
            # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd [id]
            if len(fields) < 12:
                raise ValueError(f"Malformed chain header at line {lineno}: {line.rstrip()!r}")
            blocks.append(
                {
                    "t_name": fields[2],
                    "t_size": int(fields[3]),
                    "t_strand": fields[4],
                    "q_name": fields[7],
                    "q_size": int(fields[8]),
                    "q_strand": fields[9],
                }
            )
    if not blocks:
        raise ValueError(f"No chain headers found in {path}")
    return blocks


def validate_chain_against_genomes(
    blocks: list[dict[str, Any]],
    ref_len: dict[str, int],
    query_len: dict[str, int],
) -> None:
    """Ensure chain source/target chrom names and sizes match ref/query FAI."""
    t_names = {b["t_name"] for b in blocks}
    q_names = {b["q_name"] for b in blocks}

    t_in_ref = sum(1 for n in t_names if n in ref_len)
    t_in_query = sum(1 for n in t_names if n in query_len)
    q_in_query = sum(1 for n in q_names if n in query_len)
    q_in_ref = sum(1 for n in q_names if n in ref_len)

    if t_in_query > t_in_ref and q_in_ref > q_in_query:
        raise ValueError(
            "Chain appears reversed relative to --ref_fa/--query_fa: "
            "chain tName matches query and qName matches ref. "
            "Provide a chain built with the same source->target direction."
        )

    errors: list[str] = []
    for b in blocks:
        t, ts = b["t_name"], b["t_size"]
        q, qs = b["q_name"], b["q_size"]
        if t not in ref_len:
            errors.append(f"chain source chrom {t!r} not in ref FASTA index")
        elif ref_len[t] != ts:
            errors.append(f"chain source {t!r} size {ts} != ref FAI length {ref_len[t]}")
        if q not in query_len:
            errors.append(f"chain target chrom {q!r} not in query FASTA index")
        elif query_len[q] != qs:
            errors.append(f"chain target {q!r} size {qs} != query FAI length {query_len[q]}")

    if errors:
        seen: set[str] = set()
        uniq: list[str] = []
        for e in errors:
            if e not in seen:
                seen.add(e)
                uniq.append(e)
        raise ValueError(
            "Chain does not match provided genomes:\n  - " + "\n  - ".join(uniq[:20])
        )


def _parse_scalar(value: str) -> Any:
    if value == "" or value == '""' or value == "''":
        return ""
    if (value.startswith('"') and value.endswith('"')) or (
        value.startswith("'") and value.endswith("'")
    ):
        return value[1:-1]
    if value.lower() in ("true", "false"):
        return value.lower() == "true"
    if re.fullmatch(r"-?\d+", value):
        return int(value)
    return value


def _parse_sidecar_yaml(text: str) -> dict[str, Any]:
    """Parse the small chain_meta.yml subset we emit. Fail on malformed structure."""
    if not text.strip():
        raise ValueError("chain_meta sidecar is empty")

    result: dict[str, Any] = {}
    section: str | None = None
    current_pair: dict[str, str] | None = None
    pairs: list[dict[str, str]] = []
    validation: dict[str, Any] = {}
    val_section: str | None = None
    saw_any_key = False

    def flush_pair() -> None:
        nonlocal current_pair
        if current_pair is not None:
            if "source" not in current_pair or "target" not in current_pair:
                raise ValueError("Malformed chromosome_pairs entry in chain_meta sidecar")
            pairs.append(current_pair)
            current_pair = None

    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip(" "))
        stripped = raw.strip()
        saw_any_key = True

        if indent == 0 and stripped.endswith(":") and " " not in stripped.rstrip(":"):
            flush_pair()
            section = stripped[:-1]
            val_section = None
            if section == "chromosome_pairs":
                result.setdefault("chromosome_pairs", [])
            elif section == "validation":
                result.setdefault("validation", validation)
            elif section in ("source_fasta", "target_fasta", "chain"):
                result.setdefault(section, {})
            else:
                result[section] = None
            continue

        if indent == 0 and ":" in stripped:
            flush_pair()
            section = None
            key, _, value = stripped.partition(":")
            result[key.strip()] = _parse_scalar(value.strip())
            continue

        if section == "chromosome_pairs":
            if stripped.startswith("- "):
                flush_pair()
                current_pair = {}
                rest = stripped[2:].strip()
                if rest and ":" in rest:
                    k, _, v = rest.partition(":")
                    current_pair[k.strip()] = _parse_scalar(v.strip())
                continue
            if current_pair is not None and ":" in stripped:
                k, _, v = stripped.partition(":")
                current_pair[k.strip()] = _parse_scalar(v.strip())
            continue

        if section == "validation":
            if indent == 2 and stripped.endswith(":") and not stripped.startswith("-"):
                val_section = stripped[:-1]
                validation[val_section] = {}
                continue
            if ":" in stripped:
                k, _, v = stripped.partition(":")
                key = k.strip()
                val = _parse_scalar(v.strip())
                if val_section and isinstance(validation.get(val_section), dict) and indent >= 4:
                    validation[val_section][key] = val
                else:
                    validation[key] = val
                    val_section = None
            continue

        if section in ("source_fasta", "target_fasta", "chain") and ":" in stripped:
            k, _, v = stripped.partition(":")
            result.setdefault(section, {})[k.strip()] = _parse_scalar(v.strip())
            continue

        if section is not None:
            raise ValueError(f"Malformed chain_meta sidecar line: {raw!r}")

    flush_pair()
    if not saw_any_key:
        raise ValueError("chain_meta sidecar has no recognizable keys")
    if pairs:
        result["chromosome_pairs"] = pairs
    if validation:
        result["validation"] = validation
    return result


def validate_sidecar_required(
    sidecar_path: str | None,
    ref_sha: str,
    query_sha: str,
    chain_sha: str,
) -> dict[str, Any]:
    """Require a well-formed sidecar with source/target/chain sha256 and match them."""
    if not sidecar_path:
        raise ValueError(
            "reuse with --chain_meta requires a sidecar containing "
            "source_fasta.sha256, target_fasta.sha256, and chain.sha256"
        )
    path = Path(sidecar_path)
    if not path.is_file():
        raise ValueError(f"--chain_meta file not found: {sidecar_path}")

    try:
        doc = _parse_sidecar_yaml(path.read_text(encoding="utf-8"))
    except ValueError:
        raise
    except Exception as exc:
        raise ValueError(f"Malformed chain_meta sidecar: {exc}") from exc

    for section in ("source_fasta", "target_fasta", "chain"):
        if section not in doc or not isinstance(doc[section], dict):
            raise ValueError(f"chain_meta sidecar missing required section: {section}")
        sha = doc[section].get("sha256")
        if not sha or not isinstance(sha, str):
            raise ValueError(f"chain_meta sidecar missing {section}.sha256")
        if not _SHA256_RE.fullmatch(sha):
            raise ValueError(f"chain_meta sidecar has malformed {section}.sha256: {sha!r}")

    src_sha = doc["source_fasta"]["sha256"]
    tgt_sha = doc["target_fasta"]["sha256"]
    side_chain_sha = doc["chain"]["sha256"]

    if src_sha != ref_sha:
        raise ValueError(
            f"chain_meta source_fasta.sha256 {src_sha} != current ref_fa {ref_sha}"
        )
    if tgt_sha != query_sha:
        raise ValueError(
            f"chain_meta target_fasta.sha256 {tgt_sha} != current query_fa {query_sha}"
        )
    if side_chain_sha != chain_sha:
        raise ValueError(
            f"chain_meta chain.sha256 {side_chain_sha} != provided chain {chain_sha}"
        )

    return {
        "source_sha256": src_sha,
        "target_sha256": tgt_sha,
        "chain_sha256": side_chain_sha,
        "path": sidecar_path,
    }


def ystr(value: Any) -> str:
    return json.dumps("" if value is None else str(value))


def build_chain_meta_doc(
    *,
    mode: str,
    ref_path: str,
    query_path: str,
    chain_path: str,
    ref_sha: str,
    ref_size: int,
    query_sha: str,
    query_size: int,
    chain_sha: str,
    chain_size: int,
    blocks: list[dict[str, Any]],
    pair_lines: list[str],
    pair_strategy: str,
    align_mode: str,
    split_threshold: int,
    split_size: int,
    minimap2_args: str,
    sidecar_info: dict[str, Any] | None,
) -> list[str]:
    pairs: list[dict[str, str]] = []
    for pl in pair_lines:
        parts = pl.split("\t")
        if len(parts) >= 2:
            pairs.append({"source": parts[0], "target": parts[1]})

    lines = [
        f"mode: {ystr(mode)}",
        "source_fasta:",
        f"  path: {ystr(ref_path)}",
        f"  sha256: {ystr(ref_sha)}",
        f"  bytes: {ref_size}",
        "target_fasta:",
        f"  path: {ystr(query_path)}",
        f"  sha256: {ystr(query_sha)}",
        f"  bytes: {query_size}",
        "chain:",
        f"  path: {ystr(chain_path)}",
        f"  sha256: {ystr(chain_sha)}",
        f"  bytes: {chain_size}",
        f"  blocks: {len(blocks)}",
        f"pair_strategy: {ystr(pair_strategy)}",
        f"align_mode: {ystr(align_mode)}",
        f"split_threshold: {split_threshold}",
        f"split_size: {split_size}",
        f"minimap2_args: {ystr(minimap2_args)}",
        "chromosome_pairs:",
    ]
    if pairs:
        for p in pairs:
            lines.append(f"  - source: {ystr(p['source'])}")
            lines.append(f"    target: {ystr(p['target'])}")
    else:
        lines.append("  []")
    lines.append("validation:")
    lines.append("  chain_headers_checked: true")
    lines.append(f"  orientation: {ystr('source_is_tName_target_is_qName')}")
    if sidecar_info:
        lines.append("  sidecar:")
        lines.append(f"    path: {ystr(sidecar_info.get('path'))}")
        lines.append(f"    source_sha256: {ystr(sidecar_info.get('source_sha256'))}")
        lines.append(f"    target_sha256: {ystr(sidecar_info.get('target_sha256'))}")
        if sidecar_info.get("chain_sha256"):
            lines.append(f"    chain_sha256: {ystr(sidecar_info.get('chain_sha256'))}")
    return lines


def write_chain_meta_from_files(
    *,
    ref_fa: Path,
    query_fa: Path,
    ref_fai: Path,
    query_fai: Path,
    chain: Path,
    chrom_pairs: Path,
    out_meta: Path,
    out_chain: Path,
    meta: dict[str, Any],
    require_sidecar: bool = False,
) -> None:
    mode = str(meta.get("mode") or "built")
    ref_sha, ref_size = sha256_file(ref_fa)
    query_sha, query_size = sha256_file(query_fa)
    chain_sha, chain_size = sha256_file(chain)
    ref_len = load_fai(ref_fai)
    query_len = load_fai(query_fai)
    blocks = parse_chain_headers(chain)

    validate_chain_against_genomes(blocks, ref_len, query_len)

    sidecar_info = None
    if mode == "reuse":
        sidecar = meta.get("chain_meta")
        if require_sidecar or sidecar:
            # Provided sidecar must be complete and match current inputs.
            sidecar_info = validate_sidecar_required(
                sidecar, ref_sha, query_sha, chain_sha
            )

    pair_lines: list[str] = []
    if chrom_pairs.exists() and chrom_pairs.stat().st_size > 0:
        for line in chrom_pairs.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                pair_lines.append(line)

    if mode == "reuse" and pair_lines:
        chain_pairs = {(b["t_name"], b["q_name"]) for b in blocks}
        missing = []
        for pl in pair_lines:
            parts = pl.split("\t")
            if len(parts) >= 2 and (parts[0], parts[1]) not in chain_pairs:
                missing.append(f"{parts[0]} -> {parts[1]}")
        if missing:
            raise ValueError(
                "Chromosome pairs not covered by reused chain:\n  - "
                + "\n  - ".join(missing)
            )

    minimap2_args = str(meta.get("minimap2_args") or "-cx asm5 --cs")
    lines = build_chain_meta_doc(
        mode=mode,
        ref_path=meta.get("ref_fa") or "",
        query_path=meta.get("query_fa") or "",
        chain_path=meta.get("chain_path") or "all.chain",
        ref_sha=ref_sha,
        ref_size=ref_size,
        query_sha=query_sha,
        query_size=query_size,
        chain_sha=chain_sha,
        chain_size=chain_size,
        blocks=blocks,
        pair_lines=pair_lines,
        pair_strategy=meta.get("pair_strategy") or "",
        align_mode=meta.get("align_mode") or "",
        split_threshold=int(meta.get("split_threshold") or 0),
        split_size=int(meta.get("split_size") or 0),
        minimap2_args=minimap2_args,
        sidecar_info=sidecar_info,
    )
    out_meta.write_text("\n".join(lines) + "\n", encoding="utf-8")
    # Stream copy; do not load the whole chain into memory.
    shutil.copyfile(chain, out_chain)


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-fa", type=Path, required=True)
    parser.add_argument("--query-fa", type=Path, required=True)
    parser.add_argument("--ref-fai", type=Path, required=True)
    parser.add_argument("--query-fai", type=Path, required=True)
    parser.add_argument("--chain", type=Path, required=True)
    parser.add_argument("--chrom-pairs", type=Path, required=True)
    parser.add_argument("--out-meta", type=Path, required=True)
    parser.add_argument("--out-chain", type=Path, required=True)
    parser.add_argument("--meta-json", type=str, required=True)
    parser.add_argument(
        "--require-sidecar",
        action="store_true",
        help="Require and fully validate chain_meta sidecar",
    )
    args = parser.parse_args(argv)

    meta = json.loads(args.meta_json)
    try:
        write_chain_meta_from_files(
            ref_fa=args.ref_fa,
            query_fa=args.query_fa,
            ref_fai=args.ref_fai,
            query_fai=args.query_fai,
            chain=args.chain,
            chrom_pairs=args.chrom_pairs,
            out_meta=args.out_meta,
            out_chain=args.out_chain,
            meta=meta,
            require_sidecar=args.require_sidecar,
        )
    except (ValueError, OSError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
