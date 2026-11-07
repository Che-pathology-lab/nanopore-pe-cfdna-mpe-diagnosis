"""Purpose: adapt methylKit 5 kb means to the existing NGS METH input.

Inputs: hg38_region5k.txt; discovery hyper/hypo BEDs. Output: NGS_DMR.tsv.
Run from the repository root after editing paths. No model fitting is performed.
"""
import pandas as pd

INPUT_FILE = "results/07_ngs_validation/preprocessing/methylation/hg38_region5k.txt"
HYPER_BED = "results/02_nanopore_methylation/hyper_DMR.bed"
HYPO_BED = "results/02_nanopore_methylation/hypo_DMR.bed"
OUTPUT_FILE = "results/07_ngs_validation/METH/NGS_DMR.tsv"
SAMPLE_RENAME_FILE = ""  # Optional private CSV with ID,new_name; never commit it.


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
    samples = [column for column in data if column.endswith("_meth")]
    result = data[["chr", "start", "end"] + samples].copy()
    result.columns = ["chr", "start", "end"] + [column[:-5] for column in samples]
    result["chr"] = "chr" + result["chr"].astype(str).str.removeprefix("chr")
    result["start"] = result["start"] - 1  # methylKit 1-based -> BED 0-based.
    result.iloc[:, 3:] = result.iloc[:, 3:].astype(float) * 100
    beds = [pd.read_csv(path, sep="\t", header=None).iloc[:, :3] for path in (HYPER_BED, HYPO_BED)]
    regions = pd.concat(beds).drop_duplicates()
    regions.columns = ["chr", "start", "end"]
    result = regions.merge(result, on=["chr", "start", "end"], how="inner", validate="one_to_one")
    result = rename_samples(result)
    result.to_csv(OUTPUT_FILE, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
