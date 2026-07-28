nextflow.enable.types = true

process WRITE_RUN_META {
    tag 'run_meta'
    label 'tool_py'
    label 'small_mem'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'run_meta.yml'

    input:
    versions_yml: Path
    meta_json: String

    output:
    run_meta: Path = file('run_meta.yml')

    script:
    // Base64 avoids quoting issues when embedding JSON into the task script.
    def metaB64 = meta_json.bytes.encodeBase64().toString()
    """
    python3 << 'PY'
import base64
import hashlib
import json
from pathlib import Path

meta = json.loads(base64.b64decode("${metaB64}").decode("utf-8"))
versions = Path("${versions_yml}")
sha = ""
if versions.exists() and versions.stat().st_size:
    sha = hashlib.sha256(versions.read_bytes()).hexdigest()

def dump(value):
    if value is None:
        return '""'
    return json.dumps(str(value))

lines = [
    'pipeline: "nf-liftover"',
    f"manifest_version: {dump(meta.get('manifest_version'))}",
    f"git_commit: {dump(meta.get('git_commit'))}",
    f"git_dirty: {dump(meta.get('git_dirty'))}",
    f"revision: {dump(meta.get('revision'))}",
    f"run_name: {dump(meta.get('run_name'))}",
    f"profile: {dump(meta.get('profile'))}",
    f"command_line: {dump(meta.get('command_line'))}",
    "inputs:",
    f"  id: {dump(meta.get('id'))}",
    f"  vcf: {dump(meta.get('vcf'))}",
    f"  ref_fa: {dump(meta.get('ref_fa'))}",
    f"  query_fa: {dump(meta.get('query_fa'))}",
    f"  mapping: {dump(meta.get('mapping'))}",
    f"  chain: {dump(meta.get('chain'))}",
    f"  chain_meta: {dump(meta.get('chain_meta'))}",
    "parameters:",
    f"  pair_strategy: {dump(meta.get('pair_strategy'))}",
    f"  aligner: {dump(meta.get('aligner'))}",
    f"  align_mode: {dump(meta.get('align_mode'))}",
    f"  split_threshold: {dump(meta.get('split_threshold'))}",
    f"  split_size: {dump(meta.get('split_size'))}",
    f"  flank: {dump(meta.get('flank'))}",
    f"  publish_paf: {dump(meta.get('publish_paf'))}",
    f"  outdir: {dump(meta.get('outdir'))}",
    'software_versions_file: "software_versions.yml"',
]
if sha:
    lines.append(f'software_versions_sha256: "{sha}"')

Path("run_meta.yml").write_text("\\n".join(lines) + "\\n", encoding="utf-8")
PY
    """

    stub:
    """
    echo 'pipeline: "nf-liftover"' > run_meta.yml
    """
}
