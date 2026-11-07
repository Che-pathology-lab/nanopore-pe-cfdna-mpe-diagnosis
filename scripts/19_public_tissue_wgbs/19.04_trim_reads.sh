#!/usr/bin/env bash
# Purpose: WGBS full-depth stage trim. Edit 19.00_settings.sh first.
HERE="$(cd "$(dirname "$0")" && pwd)"
bash "$HERE/_run_stage.sh" trim
