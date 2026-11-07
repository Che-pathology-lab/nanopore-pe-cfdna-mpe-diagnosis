#!/usr/bin/env bash
# Internal launcher: exports the settings and calls the preserved WGBS core.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
set -a
source "$HERE/19.00_settings.sh"
set +a
exec "$PYTHON" "$HERE/_pipeline.py" "${1:-plan}"
