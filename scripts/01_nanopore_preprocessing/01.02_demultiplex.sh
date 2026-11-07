#!/usr/bin/env bash

# Purpose: demultiplex a basecalled BAM using the original kit setting.
# Edit the example paths below before running: bash 01.02_demultiplex.sh

INPUT_BAM="results/01_nanopore_preprocessing/1.1_basecall/temp_calls.bam"
OUTPUT_DIR="results/01_nanopore_preprocessing/1.2_demux"
KIT_NAME="SQK-NBD114-96"
THREADS=16

dorado demux \
    --no-classify \
    -t "$THREADS" \
    --kit-name "$KIT_NAME" \
    --output-dir "$OUTPUT_DIR" \
    "$INPUT_BAM"
