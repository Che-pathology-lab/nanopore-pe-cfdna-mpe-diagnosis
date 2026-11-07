#!/usr/bin/env bash

# Purpose: generate 5 Mb readCounter WIG files for ichorCNA.
# Edit the example paths below before running: bash 03.01_wigmake.sh

INPUT_DIR="results/01_nanopore_preprocessing/1.4_mapping"
OUTPUT_DIR="results/03_nanopore_cna/3.1_wig"
MAPQ=20
CHROMOSOMES="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

for bam_file in "$INPUT_DIR"/*.bam; do
    base_name="$(basename "$bam_file" .bam)"
    output_wig="${OUTPUT_DIR}/${base_name}.wig"
    readCounter --window 5000000 --quality "$MAPQ" \
        --chromosome "$CHROMOSOMES" \
        "$bam_file" > "$output_wig"
done
