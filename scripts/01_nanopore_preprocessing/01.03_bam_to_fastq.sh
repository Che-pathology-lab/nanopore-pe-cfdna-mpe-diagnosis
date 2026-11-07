#!/usr/bin/env bash

# Purpose: convert one demultiplexed barcode BAM to a quality-filtered FASTQ.
# Run this script separately for each barcode BAM.

INPUT_BAM="results/01_nanopore_preprocessing/1.2_demux/barcode01.bam"
OUTPUT_DIR="results/01_nanopore_preprocessing/1.3_fastq"

# Set SAMPLE_NAME to the corresponding pseudonymized sample name according to
# the barcode-to-sample mapping table.
SAMPLE_NAME="sample_name"

THREADS=10

samtools fastq -T MM,ML -@ "$THREADS" "$INPUT_BAM" | \
    NanoFilt -q 7 | \
    pigz --best > "${OUTPUT_DIR}/${SAMPLE_NAME}.fastq.gz"