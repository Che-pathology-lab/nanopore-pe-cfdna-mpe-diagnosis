#!/usr/bin/env bash

# Purpose: summarize methylation calls in the original 5 kb regions.
# Edit the example paths below before running: bash 02.02_5kb_bin_modstat.sh

INPUT_DIR="results/02_nanopore_methylation/2.1_modkit_file"
OUTPUT_DIR="results/02_nanopore_methylation/2.2_modkit_stats"
BED_FILE="references/hg38_genome_wide_5kb_windows.bed"
THREADS=6

for bed_file in "$INPUT_DIR"/*.bed.gz; do
    base_name="$(basename "$bed_file" .bed.gz)"
    output_file="${OUTPUT_DIR}/${base_name}.tsv"
    modkit stats -t "$THREADS" "$bed_file" \
        --regions "$BED_FILE" \
        -o "$output_file"
done
