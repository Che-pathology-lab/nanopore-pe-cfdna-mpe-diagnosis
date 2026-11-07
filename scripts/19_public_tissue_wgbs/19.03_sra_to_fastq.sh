#!/usr/bin/env bash
# Purpose: WGBS full-depth stage fastq. Edit 19.00_settings.sh first.
HERE="$(cd "$(dirname "$0")" && pwd)"
bash "$HERE/_run_stage.sh" fastq
