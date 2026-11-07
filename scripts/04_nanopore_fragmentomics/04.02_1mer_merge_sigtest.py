"""Purpose: 1mer merge sigtest for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import os
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


INPUT_FILE = "results/04_nanopore_fragmentomics/4mer_endmotif_matrix.csv"
GROUP_PATH = "metadata/discovery_groupfile.csv"
CSV_OUTPUT_NAME = "results/04_nanopore_fragmentomics/motif_1mer_significant.csv"

FDR_THRESHOLD = 0.01


def main():
    os.makedirs(os.path.dirname(CSV_OUTPUT_NAME), exist_ok=True)

    df_features = pd.read_csv(INPUT_FILE, index_col=0)
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

    motif_cols = [
        c for c in df_features.columns
        if len(str(c)) == 4 and set(str(c)).issubset(set("ACGT"))
    ]

    onemer_features = ["A", "C", "G", "T"]

    df_1mer = merged_df[["Sample", "group", "split"]].copy()

    for monomer in onemer_features:
        cols_to_sum = [
            c for c in motif_cols
            if str(c).startswith(monomer)
        ]

        if cols_to_sum:
            df_1mer[monomer] = merged_df[cols_to_sum].apply(
                pd.to_numeric,
                errors="coerce",
            ).sum(axis=1, skipna=True)

    train_df = df_1mer[
        (df_1mer["split"] == "train")
        & (df_1mer["group"].isin(["mpe", "bpe"]))
    ].copy()

    stats_results = []

    for feature in onemer_features:
        vals_mpe = train_df.loc[
            train_df["group"] == "mpe",
            feature,
        ].dropna()

        vals_bpe = train_df.loc[
            train_df["group"] == "bpe",
            feature,
        ].dropna()

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

    _, pvals_corrected, _, _ = multipletests(
        res_df["P_Raw"],
        alpha=FDR_THRESHOLD,
        method="fdr_bh",
    )

    res_df["P_Adj"] = pvals_corrected
    res_df["Significant"] = res_df["P_Adj"] < FDR_THRESHOLD

    sig_features = res_df.loc[
        res_df["Significant"] == True,
        "Feature",
    ].tolist()

    if not sig_features:
        df_1mer[["Sample"]].to_csv(CSV_OUTPUT_NAME, index=False)
    else:
        cols_to_save = ["Sample"] + sig_features
        df_1mer[cols_to_save].to_csv(CSV_OUTPUT_NAME, index=False)



if __name__ == "__main__":
    main()
