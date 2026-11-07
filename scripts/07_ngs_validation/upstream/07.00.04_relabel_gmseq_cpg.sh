#!/usr/bin/env bash

# Purpose: relabel GM-seq converted-C fractions and shift unmatched reverse sites.
# Create OUTPUT_DIR before running. The merged files are not edited in place.
INPUT_DIR="results/07_ngs_validation/preprocessing/cpg_merged"
OUTPUT_DIR="results/07_ngs_validation/preprocessing/cpg_methylkit"

for input_file in "$INPUT_DIR"/*.CpG.txt; do
    filename=$(basename "$input_file")
    awk '
        NR == 1 { $6 = "freqC"; $7 = "freqT" }
        NR > 1 && $4 == "R" { $3 = $3 - 1; $4 = "F" }
        { print }
    ' OFS="\t" "$input_file" > "${OUTPUT_DIR}/${filename}"
done
