"""Purpose: merge conumee shift ratios into the existing CHIP CNA input.

Inputs: per-sample .bins.txt files, optionally a private ID/new_name mapping.
Output: CHIP_shift_ratio.csv; coordinates and ratios are retained unchanged.
Run from the repository root after editing paths.
"""
import os

import pandas as pd

INPUT_DIR = "results/08_chip_validation/preprocessing/cna/bins"
OUTPUT_FILE = "results/08_chip_validation/CHIP_CNA/CHIP_shift_ratio.csv"
SAMPLE_RENAME_FILE = ""  # Leave empty when array metadata already uses public names.
BIN_TEMPLATE_FILE = "references/chip_cna_bin_template.tsv"


def main():
    sample_map = None
    if SAMPLE_RENAME_FILE:
        mapping = pd.read_csv(SAMPLE_RENAME_FILE, dtype=str)
        if mapping["ID"].duplicated().any() or mapping["new_name"].duplicated().any():
            raise ValueError("Sample mapping must contain unique ID and new_name values.")
        sample_map = mapping.set_index("ID")["new_name"].to_dict()
    result = None
    coordinates = ["chr", "start", "end", "probes"]
    for filename in sorted(os.listdir(INPUT_DIR)):
        if not filename.endswith(".bins.txt"):
            continue
        sample = filename[:-9]
        if sample_map is not None:
            sample = sample_map[sample]
        data = pd.read_csv(os.path.join(INPUT_DIR, filename), sep="\t")
        data = data.rename(columns={"seqnames": "chr", "shift_ratio": sample})
        data = data[coordinates + [sample]]
        if result is None:
            result = data
        else:
            if sample in result.columns:
                raise ValueError("Duplicate sample after mapping.")
            result = result.merge(data, on=coordinates, how="outer", validate="one_to_one")
    if result is None:
        raise ValueError("No .bins.txt input files found.")
    # Freeze the published 576-bin coordinate set. The supplied conumee tables
    # contain 616 bins; their historical 40-bin exclusion rule was not supplied.
    # An empty template path exports all native conumee bins for inspection.
    if BIN_TEMPLATE_FILE:
        template = pd.read_csv(BIN_TEMPLATE_FILE, sep="\t")
        result = template.merge(result, on=coordinates, how="left", validate="one_to_one")
    result.to_csv(OUTPUT_FILE, index=False)


if __name__ == "__main__":
    main()
