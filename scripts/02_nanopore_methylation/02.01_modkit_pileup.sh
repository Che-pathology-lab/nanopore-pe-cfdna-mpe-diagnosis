#!/usr/bin/env bash

# Purpose: generate traditional-preset modkit pileups for mapped BAM files.
# Edit the example paths below before running: bash 02.01_modkit_pileup.sh

INPUT_DIR="results/01_nanopore_preprocessing/1.4_mapping"
OUTPUT_DIR="results/02_nanopore_methylation/2.1_modkit_file"
REFERENCE_FASTA="references/hg38.fa"
THREADS=6

for bam_file in "$INPUT_DIR"/*.bam; do
    base_name="$(basename "$bam_file" .bam)"
    output_file="${OUTPUT_DIR}/${base_name}.bed.gz"
    modkit pileup "$bam_file" /dev/stdout \
        --preset traditional \
        --ref "$REFERENCE_FASTA" \
        --filter-threshold C:0.7 -t "$THREADS" | \
        bgzip -l 9 -@ "$THREADS" > "$output_file"
done
