nextflow.enable.types = true

process WRITE_CHAIN_META {
    tag "${mode}"
    label 'tool_py'
    label 'small_mem'

    // Always publish metadata and a copy of the chain (including reuse runs).
    input:
    ref_fa: Path
    query_fa: Path
    ref_fai: Path
    query_fai: Path
    chrom_pairs: Path
    chain: Path
    mode: String
    meta_json: String
    record(
        chain_meta_file: Path?
    )

    stage:
    stageAs ref_fa, 'meta_ref.fa'
    stageAs query_fa, 'meta_query.fa'
    stageAs ref_fai, 'meta_ref.fa.fai'
    stageAs query_fai, 'meta_query.fa.fai'
    stageAs chrom_pairs, 'chrom_pairs.tsv'
    stageAs chain, 'input.chain'
    stageAs chain_meta_file, 'sidecar_chain_meta.yml'

    output:
    meta: Path = file('chain_meta.yml')
    chain_out: Path = file('all.chain')
    versions: Path = file('versions.yml')

    script:
    // Base64 avoids shell expansion / quoting issues for paths and JSON.
    def metaB64 = meta_json.bytes.encodeBase64().toString()
    def hasSidecar = chain_meta_file ? 'true' : 'false'
    """
    python3 - <<'PY'
import base64
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, "${projectDir}/bin")
from chain_meta import write_chain_meta_from_files

meta = json.loads(base64.b64decode("${metaB64}").decode("utf-8"))
# Prefer staged sidecar so tasks do not depend on absolute host paths.
if "${hasSidecar}" == "true" and Path("sidecar_chain_meta.yml").is_file():
    meta["chain_meta"] = "sidecar_chain_meta.yml"
try:
    write_chain_meta_from_files(
        ref_fa=Path("meta_ref.fa"),
        query_fa=Path("meta_query.fa"),
        ref_fai=Path("meta_ref.fa.fai"),
        query_fai=Path("meta_query.fa.fai"),
        chain=Path("input.chain"),
        chrom_pairs=Path("chrom_pairs.tsv"),
        out_meta=Path("chain_meta.yml"),
        out_chain=Path("all.chain"),
        meta=meta,
    )
except Exception as exc:
    print(str(exc), file=sys.stderr)
    raise SystemExit(1)
PY

    {
      echo '"${task.process}":'
      echo "    python: \$(python3 --version 2>&1 | awk '{ print \$2 }')"
      echo "    python_path: \$(command -v python3)"
    } > versions.yml
    """

    stub:
    """
    cp -L input.chain all.chain 2>/dev/null || touch all.chain
    cat > chain_meta.yml <<'EOF'
mode: "stub"
source_fasta:
  path: "stub"
  sha256: "0"
  bytes: 0
target_fasta:
  path: "stub"
  sha256: "0"
  bytes: 0
chain:
  path: "stub"
  sha256: "0"
  bytes: 0
  blocks: 0
chromosome_pairs: []
validation:
  chain_headers_checked: false
  orientation: "source_is_tName_target_is_qName"
EOF
    {
      echo '"${task.process}":'
      echo "    python: 3.11.0"
      echo "    python_path: stub"
    } > versions.yml
    """
}
