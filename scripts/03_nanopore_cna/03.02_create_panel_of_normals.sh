#!/usr/bin/env bash

# Purpose: create the ichorCNA panel of normals with the original references.
# Edit the example paths below before running: bash 03.02_create_panel_of_normals.sh

PON_LIST="references/ichorCNA/panelofnormalfile.txt"
REF_DIR="references/ichorCNA"
OUTPUT_PREFIX="results/03_nanopore_cna/panelofnormals"

Rscript createPanelOfNormals.R \
    --filelist "$PON_LIST" \
    --gcWig "${REF_DIR}/gc_hg38_5000kb_clean.wig" \
    --mapWig "${REF_DIR}/map_hg38_5000kb_clean.wig" \
    --centromere "${REF_DIR}/GRCh38.GCA_000001405.2_centromere_acen.txt" \
    --outfile "$OUTPUT_PREFIX"
