#!/usr/bin/env bash
# Run full-depth methylation/CNA only. Downloads are a separate optional step.
HERE="$(cd "$(dirname "$0")" && pwd)"
for stage in check validate fastq trim bsmap bam mcall dmr merge_dmr wig cna merge_cna qc; do
    bash "$HERE/_run_stage.sh" "$stage" || exit 1
done
