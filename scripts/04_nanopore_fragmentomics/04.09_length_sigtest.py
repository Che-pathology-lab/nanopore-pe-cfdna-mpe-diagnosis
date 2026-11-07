"""Purpose: Length sigtest for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


FEATURE_PATH = "results/04_nanopore_fragmentomics/nanopore_lengths_matrix.csv"
GROUP_PATH = "metadata/discovery_groupfile.csv"
OUTPUT_FILE = "results/04_nanopore_fragmentomics/discovery_size_feature_significant.csv"

FDR_THRESHOLD = 0.01

FRAGMENT_RANGES = [
    ("short_mono", 100, 150, "100-150 bp\nShort"),
    ("canonical_mono", 150, 170, "150-170 bp\nCanonical mono"),
    ("long_mono", 170, 220, "170-220 bp\nLong mono"),
    ("mono_di_transition", 220, 270, "220-270 bp\nMono-di transition"),
    ("di_nucleosomal", 270, 320, "Di-nucleosomal"),
]


def get_length_columns(feature_cols, start, end):
    selected_cols = []

    for col in feature_cols:
        try:
            length = int(float(str(col)))
        except ValueError:
            continue

        if start <= length <= end:
            selected_cols.append(col)

    return selected_cols


def main():
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

    df_features = pd.read_csv(FEATURE_PATH, index_col=0)
    df_labels = pd.read_csv(GROUP_PATH)

    df_features.index = df_features.index.astype(str).str.strip()
    df_labels["Sample"] = df_labels["Sample"].astype(str).str.strip()
    df_labels["group"] = df_labels["group"].astype(str).str.strip().str.lower()
    df_labels["split"] = df_labels["split"].astype(str).str.strip().str.lower()

    df_labels["split"] = df_labels["split"].replace({
        "panelofnormal": "train",
    })

    merged_df = df_labels.merge(
        df_features,
        left_on="Sample",
        right_index=True,
        how="inner",
    )

    feature_cols = df_features.columns.tolist()
    size_feature_names = []

    for feature_id, start, end, label in FRAGMENT_RANGES:
        target_cols = get_length_columns(
            feature_cols=feature_cols,
            start=start,
            end=end,
        )

        if not target_cols:
            continue

        feature_name = f"Freq_{start}_{end}"
        size_feature_names.append(feature_name)

        merged_df[feature_name] = merged_df[target_cols].apply(
            pd.to_numeric,
            errors="coerce",
        ).sum(axis=1, skipna=True)

    train_df = merged_df[
        (merged_df["split"] == "train")
        & (merged_df["group"].isin(["mpe", "bpe"]))
    ].copy()

    stats_results = []

    for feature in size_feature_names:
        vals_mpe = train_df.loc[
            train_df["group"] == "mpe",
            feature,
        ].dropna()

        vals_bpe = train_df.loc[
            train_df["group"] == "bpe",
            feature,
        ].dropna()

        if len(vals_mpe) == 0 or len(vals_bpe) == 0:
            p_raw = np.nan
        else:
            _, p_raw = mannwhitneyu(
                vals_mpe,
                vals_bpe,
                alternative="two-sided",
            )

        stats_results.append({
            "Feature": feature,
            "P_Raw": p_raw,
        })

    res_df = pd.DataFrame(stats_results)
    res_df["P_Adj"] = np.nan
    res_df["Significant"] = False

    valid_mask = res_df["P_Raw"].notna()

    if valid_mask.any():
        _, pvals_corrected, _, _ = multipletests(
            res_df.loc[valid_mask, "P_Raw"],
            alpha=FDR_THRESHOLD,
            method="fdr_bh",
        )

        res_df.loc[valid_mask, "P_Adj"] = pvals_corrected
        res_df.loc[valid_mask, "Significant"] = (
            res_df.loc[valid_mask, "P_Adj"] < FDR_THRESHOLD
        )

    sig_features = res_df.loc[
        res_df["Significant"] == True,
        "Feature",
    ].tolist()

    if not sig_features:
        merged_df[["Sample"]].to_csv(OUTPUT_FILE, index=False)
    else:
        merged_df[["Sample"] + sig_features].to_csv(OUTPUT_FILE, index=False)


if __name__ == "__main__":
    main()