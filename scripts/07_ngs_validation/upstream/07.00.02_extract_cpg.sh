#!/usr/bin/env bash

# Purpose: extract GM-seq cytosine reports and prepare seven-column CpG tables.
# Create each chromosome directory and OUTPUT_DIR before running.
SAMPLE_NAME="sample_name"
INPUT_BAM="results/07_ngs_validation/preprocessing/bam/sample_name_sorted.bam"
REFERENCE_FASTA="references/GRCh38.fa"
CHROM_LENGTHS="references/GRCh38.chrom.sizes"
REPORT_DIR="results/07_ngs_validation/preprocessing/cytosine_reports"
OUTPUT_DIR="results/07_ngs_validation/preprocessing/cpg_raw"
THREADS=1
CHROMS=(chr{1..22} chrX chrY)

# Confirm this suffix against the installed MethylDackel build.
CYTOSINE_SUFFIX=".cytosine_report.txt"
# gctk must be available on PATH; its source/version was not supplied.
printf 'chrBase\tchr\tbase\tstrand\tcoverage\tfreqT\tfreqC\n' \
    > "${OUTPUT_DIR}/${SAMPLE_NAME}.CpG.txt"

for chrom in "${CHROMS[@]}"; do
    report_prefix="${REPORT_DIR}/${chrom}/${SAMPLE_NAME}_${chrom}_cpg"
    cpg_file="${REPORT_DIR}/${chrom}/${SAMPLE_NAME}_${chrom}_CpG.txt"
    MethylDackel extract --mCtoT --cytosine_report --CHG --CHH \
        "$REFERENCE_FASTA" "$INPUT_BAM" -p 20 -q 10 -@ "$THREADS" \
        -r "$chrom" -o "$report_prefix"
    gctk prepareCpG -s "$SAMPLE_NAME" -c "$chrom" \
        -m "${report_prefix}${CYTOSINE_SUFFIX}" \
        -l "$CHROM_LENGTHS" -o "$cpg_file"
    # Retain the original concatenation behavior. Confirm that gctk
    # writes data rows only, without a repeated chromosome header.
    cat "$cpg_file" >> "${OUTPUT_DIR}/${SAMPLE_NAME}.CpG.txt"
done
