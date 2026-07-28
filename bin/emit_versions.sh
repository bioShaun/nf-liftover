#!/usr/bin/env bash
# Emit a minimal valid YAML versions fragment without heredoc delimiter issues.
# Usage: emit_versions.sh "PROCESS_NAME" key=value [key=value ...]
set -euo pipefail
process_name="${1:?process name required}"
shift
echo "\"${process_name}\":"
if [ "$#" -eq 0 ]; then
  echo "    status: ok"
  exit 0
fi
for kv in "$@"; do
  key="${kv%%=*}"
  val="${kv#*=}"
  printf '    %s: %s\n' "${key}" "${val}"
done
