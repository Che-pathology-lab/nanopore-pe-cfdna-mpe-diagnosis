"""Purpose: adapt EPIC v2 DMR beta means to the existing CHIP METH input.

Inputs: chip_dmr_beta.tsv and an approved private ID/new_name mapping.
Output: CHIP_DMR.tsv with 0-based coordinates and methylation percentages.
Run from the repository root after editing paths.
"""
import pandas as pd

INPUT_FILE = "results/08_chip_validation/preprocessing/methylation/chip_dmr_beta.tsv"
OUTPUT_FILE = "results/08_chip_validation/CHIP_METH/CHIP_DMR.tsv"
SAMPLE_RENAME_FILE = "metadata/chip_preprocessing/sample_rename.csv"


def rename_samples(frame):
    if SAMPLE_RENAME_FILE:
        mapping = pd.read_csv(SAMPLE_RENAME_FILE, dtype=str)
        if mapping["ID"].duplicated().any() or mapping["new_name"].duplicated().any():
            raise ValueError("Sample mapping must contain unique ID and new_name values.")
        names = mapping.set_index("ID")["new_name"].to_dict()
        missing = [sample for sample in frame.columns[3:] if sample not in names]
        if missing:
            raise ValueError("Some source samples lack a private ID/new_name mapping.")
        frame = frame.rename(columns=names)
    if frame.columns.duplicated().any():
        raise ValueError("Duplicate sample columns after conversion.")
    return frame


def main():
    data = pd.read_csv(INPUT_FILE, sep="\t")
    samples = [column for column in data if column not in ["ID", "chr", "start", "end"]]
    result = data[["chr", "start", "end"] + samples].copy()
    result.iloc[:, 3:] = result.iloc[:, 3:].astype(float) * 100
    result = rename_samples(result)
    result.to_csv(OUTPUT_FILE, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
