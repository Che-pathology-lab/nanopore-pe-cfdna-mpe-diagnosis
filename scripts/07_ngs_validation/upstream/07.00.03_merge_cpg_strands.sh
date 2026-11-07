#!/usr/bin/env bash

# Purpose: combine forward/reverse GM-seq CpG counts using the supplied AWK rules.
# Create OUTPUT_DIR before running. Reads INPUT_DIR/*.CpG.txt.
INPUT_DIR="results/07_ngs_validation/preprocessing/cpg_raw"
OUTPUT_DIR="results/07_ngs_validation/preprocessing/cpg_merged"

# Original strand matching and rounding are retained, including row-order behavior.
for input_file in "$INPUT_DIR"/*.CpG.txt; do
    filename=$(basename "$input_file")
    output_file="${OUTPUT_DIR}/${filename}"
    head -n 1 "$input_file" > "$output_file"
    awk -F'\t' -v OFS='\t' -v output_file="$output_file" '
    BEGIN {

        delete processed
    }


    NR == FNR {
        if (NR == 1) next;

        col1 = $1
        chr = $2
        base = $3
        strand = $4
        cov = $5
        freqT = $6
        freqC = $7

        if (strand == "R") {

            key = chr ":" base

            r_col1[key] = col1
            r_cov[key] = cov
            r_freqT[key] = freqT
            r_freqC[key] = freqC
        }
        next
    }


    {
        if (FNR == 1) next;

        col1 = $1
        chr = $2
        base = $3
        strand = $4
        cov = $5
        freqT = $6
        freqC = $7


        if (processed[FNR]) next

        if (strand == "F") {

            target_base = base + 1
            key = chr ":" target_base


            if (key in r_col1) {

                merged_cov = cov + r_cov[key]


                t1 = int(freqT * 0.01 * cov + 0.5)
                t2 = int(r_freqT[key] * 0.01 * r_cov[key] + 0.5)
                sum_t = t1 + t2
                merged_freqT = (sum_t / merged_cov) * 100


                c1 = int(freqC * 0.01 * cov + 0.5)
                c2 = int(r_freqC[key] * 0.01 * r_cov[key] + 0.5)
                sum_c = c1 + c2
                merged_freqC = (sum_c / merged_cov) * 100


                printf("%s\t%s\t%s\t%s\t%d\t%.2f\t%.2f\n",
                       col1, chr, base, strand,
                       merged_cov, merged_freqT, merged_freqC) >> output_file


                processed[key] = 1
            } else {

                print col1, chr, base, strand, cov, freqT, freqC >> output_file
            }
        } else if (strand == "R") {

            key = chr ":" base
            if (!(key in processed)) {

                split(col1, parts, /\./)
                if (length(parts) == 2) {

                    num = parts[2] + 0
                    new_num = num - 1
                    new_col1 = parts[1] "." new_num
                } else {

                    new_col1 = col1
                }

                print new_col1, chr, base, strand, cov, freqT, freqC >> output_file
            }
        } else {

            print col1, chr, base, strand, cov, freqT, freqC >> output_file
        }
    }' "$input_file" "$input_file"


done
