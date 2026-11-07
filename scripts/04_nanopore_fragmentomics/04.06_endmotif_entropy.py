"""Purpose: Endmotif entropy for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import entropy, mannwhitneyu

INPUT_FILE = "results/04_nanopore_fragmentomics/4mer_endmotif_matrix.csv"
GROUP_PATH = "metadata/discovery_groupfile.csv"

OUTPUT_CSV = "results/04_nanopore_fragmentomics/motif_entropy_matrix.csv"
OUTPUT_STATS_CSV = "results/04_nanopore_fragmentomics/motif_entropy_train_mannwhitneyu.csv"


def calculate_entropy():
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)

    df_motif = pd.read_csv(INPUT_FILE, index_col=0)
    df_motif.index = df_motif.index.astype(str).str.strip()

    df_group = pd.read_csv(GROUP_PATH)
    df_group["Sample"] = df_group["Sample"].astype(str).str.strip()
    df_group["group"] = df_group["group"].astype(str).str.strip().str.lower()
    df_group["split"] = df_group["split"].astype(str).str.strip().str.lower()
    df_group["split"] = df_group["split"].replace({
        "panelofnormal": "train",
    })

    motif_cols = [
        c for c in df_motif.columns
        if len(str(c)) == 4 and set(str(c)).issubset({"A", "C", "G", "T"})
    ]

    entropy_list = []

    for _, row in df_motif.iterrows():
        counts = row[motif_cols].values.astype(float)
        total = counts.sum()

        if total <= 0:
            ent = np.nan
        else:
            probs = counts / total
            ent = entropy(probs, base=2)

        entropy_list.append(ent)

    df_entropy = df_motif[[]].copy()
    df_entropy["Motif_Entropy"] = entropy_list
    df_entropy.reset_index(inplace=True)
    df_entropy.columns = ["Sample", "Motif_Entropy"]

    df_entropy = df_group[["Sample", "group", "split"]].merge(
        df_entropy,
        on="Sample",
        how="inner",
    )

    df_entropy[["Sample", "Motif_Entropy"]].to_csv(
        OUTPUT_CSV,
        index=False,
    )

    train_df = df_entropy[
        (df_entropy["split"] == "train")
        & (df_entropy["group"].isin(["mpe", "bpe"]))
    ].copy()

    vals_mpe = train_df.loc[
        train_df["group"] == "mpe",
        "Motif_Entropy",
    ].dropna()

    vals_bpe = train_df.loc[
        train_df["group"] == "bpe",
        "Motif_Entropy",
    ].dropna()

    if len(vals_mpe) == 0 or len(vals_bpe) == 0:
        stat = np.nan
        p_value = np.nan
    else:
        stat, p_value = mannwhitneyu(
            vals_mpe,
            vals_bpe,
            alternative="two-sided",
        )

    stats_df = pd.DataFrame([{
        "Feature": "Motif_Entropy",
        "mpe_median": vals_mpe.median(),
        "bpe_median": vals_bpe.median(),
        "p_value": p_value,
    }])

    stats_df.to_csv(
        OUTPUT_STATS_CSV,
        index=False,
    )


def main() -> None:
    calculate_entropy()


if __name__ == "__main__":
    main()
