#!/usr/bin/env bash
# Optional download template. Skip when SRA files are already available.
# No downloads are performed merely by integrating this module.
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/19.00_settings.sh"
RUN_LIST="metadata/public_tissue_wgbs/run_accessions.txt"
while IFS= read -r run_id; do
    "${SRA_BIN}/prefetch" "$run_id" --output-directory "$SRA_DIR" || exit 1
done < "$RUN_LIST"
