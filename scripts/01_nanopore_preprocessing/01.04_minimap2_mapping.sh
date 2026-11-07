#!/usr/bin/env bash

# Purpose: map Nanopore FASTQ files and create sorted, indexed BAM files.

INPUT_DIR="results/01_nanopore_preprocessing/1.3_fastq"
REFERENCE_MMI="references/hg38.mmi"
OUTPUT_DIR="results/01_nanopore_preprocessing/1.4_mapping"
MINIMAP_THREADS=8
SORT_THREADS=8

for fastq_file in "$INPUT_DIR"/*.fastq.gz; do
    sample_name=$(basename "$fastq_file" .fastq.gz)
    final_bam="${OUTPUT_DIR}/${sample_name}.bam"

    minimap2 -x map-ont -a -t "$MINIMAP_THREADS" -y --secondary=no \
        "$REFERENCE_MMI" "$fastq_file" | \
        samtools sort -@ "$SORT_THREADS" -O BAM -o "$final_bam"

    samtools index "$final_bam"
done