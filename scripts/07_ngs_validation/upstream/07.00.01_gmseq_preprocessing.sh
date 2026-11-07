#!/usr/bin/env bash

# Purpose: trim paired GM-seq reads and create a coordinate-sorted BAM.
# Run from the repository root. Create output directories before running.
SAMPLE_NAME="sample_name"
RAW_DIR="data/ngs/raw"
CLEAN_DIR="results/07_ngs_validation/preprocessing/clean_reads"
BAM_DIR="results/07_ngs_validation/preprocessing/bam"
REFERENCE_FASTA="references/GRCh38.fa"
THREADS=20

fastp \
    -i "${RAW_DIR}/${SAMPLE_NAME}_raw_1.fq.gz" \
    -I "${RAW_DIR}/${SAMPLE_NAME}_raw_2.fq.gz" \
    -o "${CLEAN_DIR}/${SAMPLE_NAME}_clean_1.fq.gz" \
    -O "${CLEAN_DIR}/${SAMPLE_NAME}_clean_2.fq.gz" \
    -h "${CLEAN_DIR}/${SAMPLE_NAME}_fastp.html" \
    -j "${CLEAN_DIR}/${SAMPLE_NAME}_fastp.json" \
    -f 6 -F 6 --length_required 35 -q 15 -u 40 -n 5

sentieon bwa mem "$REFERENCE_FASTA" \
    "${CLEAN_DIR}/${SAMPLE_NAME}_clean_1.fq.gz" \
    "${CLEAN_DIR}/${SAMPLE_NAME}_clean_2.fq.gz" | \
    sentieon util sort -o "${BAM_DIR}/${SAMPLE_NAME}_sorted.bam" \
        -t "$THREADS" --sam2bam

# The supplied code has no duplicate-marking or SNP-exclusion commands.
# Confirm these provider-side steps before treating this as a full pipeline.
samtools index "${BAM_DIR}/${SAMPLE_NAME}_sorted.bam"
