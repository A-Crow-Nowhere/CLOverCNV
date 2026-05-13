#!/usr/bin/env bash
set -euo pipefail
MAPI_ROOT="${MAPI_ROOT:-$HOME/MalariAPI}"
PY="$MAPI_ROOT/bin/modules/.CLOverCNV/cnv_post_fuse.py"
ENV="$MAPI_ROOT/envs/CLOverCNV"

if [[ -d "$ENV" ]]; then
  exec "$MAPI_ROOT/tools/miniconda3/bin/conda" run -p "$ENV" --no-capture-output python "$PY" "$@"
else
  exec python "$PY" "$@"
fi
