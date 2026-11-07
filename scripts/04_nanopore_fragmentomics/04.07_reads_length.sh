#!/usr/bin/env bash

# Purpose: calculate aligned read-length summaries with NanoPlot.
# Edit the example paths below before running: bash 04.07_reads_length.sh

INPUT_DIR="results/01_nanopore_preprocessing/1.4_mapping"
OUTPUT_DIR="results/04_nanopore_fragmentomics/nanoplot"
THREADS=4

for bam_file in "$INPUT_DIR"/*.bam; do
    base_name="$(basename "$bam_file" .bam)"
    NanoPlot \
        -t "$THREADS" \
        --bam "$bam_file" \
        --alength \
        --raw \
        --tsv_stats \
        -o "$OUTPUT_DIR" \
        -p "${base_name}_" \
        --no_static
done
