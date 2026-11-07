#!/usr/bin/env bash

# Purpose: basecall POD5 reads while retaining the original model, kit, and quality settings.
# Edit the example paths below before running: bash 01.01_basecall.sh

POD5_DIR="data/raw/pod5"
OUTPUT_DIR="results/01_nanopore_preprocessing/1.1_basecall"
MODEL_DIR="references/dorado_models"
KIT_NAME="SQK-NBD114-96"

dorado basecaller hac,5mC_5hmC \
    "$POD5_DIR" \
    --recursive \
    --models-directory "$MODEL_DIR" \
    --disable-read-splitting -x auto \
    --kit-name "$KIT_NAME" \
    --no-trim \
    --min-qscore 7 > "${OUTPUT_DIR}/temp_calls.bam"
