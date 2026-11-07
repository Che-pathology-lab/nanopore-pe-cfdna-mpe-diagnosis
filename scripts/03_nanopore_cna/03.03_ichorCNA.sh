#!/usr/bin/env bash

# Purpose: run ichorCNA with the manuscript's original model settings.
# Edit the example paths below before running: bash 03.03_ichorCNA.sh

INPUT_DIR="results/03_nanopore_cna/3.1_wig"
OUTPUT_DIR="results/03_nanopore_cna/3.3_ichorCNA"
REF_DIR="references/ichorCNA"
PON_FILE="references/ichorCNA/panelofnormals_median.rds"

for wig_file in "$INPUT_DIR"/*.wig; do
    base_name="$(basename "$wig_file" .wig)"
    Rscript runIchorCNA.R \
        --id "$base_name" \
        --WIG "$wig_file" \
        --ploidy "c(2)" \
        --normal "c(0.4,0.5,0.6,0.7,0.8,0.9)" \
        --maxCN 3 \
        --gcWig "${REF_DIR}/gc_hg38_5000kb_clean.wig" \
        --mapWig "${REF_DIR}/map_hg38_5000kb_clean.wig" \
        --centromere "${REF_DIR}/GRCh38.GCA_000001405.2_centromere_acen.txt" \
        --normalPanel "$PON_FILE" \
        --includeHOMD FALSE \
        --chrs "c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)" \
        --chrTrain "c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)" \
        --genomeStyle UCSC \
        --estimateNormal TRUE \
        --estimatePloidy TRUE \
        --estimateScPrevalence FALSE \
        --altFracThreshold 0.01 \
        --rmCentromereFlankLength 1000000 \
        --outDir "$OUTPUT_DIR"
done
