#!/usr/bin/env bash

# Purpose: compute chromosome coverage for all BAM files and generate a MultiQC report.

INPUT_DIR="results/01_nanopore_preprocessing/1.4_mapping"
BED_FILE="references/chromosomes1-22xy.bed"
OUTPUT_DIR="results/01_nanopore_preprocessing/1.5_coverage"

for bam_file in "$INPUT_DIR"/*.bam; do
    sample_name=$(basename "$bam_file" .bam)
    output_prefix="${OUTPUT_DIR}/${sample_name}"

    mosdepth -n -b "$BED_FILE" -x "$output_prefix" "$bam_file"
done

multiqc -f -o "$OUTPUT_DIR" "$OUTPUT_DIR"