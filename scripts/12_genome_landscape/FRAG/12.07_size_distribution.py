"""Purpose: Size distribution for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

LENGTH_MATRIX = "results/12_genome_landscape/FRAG/ONT_lengths_matrix.csv"
GROUP_FILE = "metadata/discovery_groupfile.csv"

OUT_DIR = "results/12_genome_landscape/FRAG"
OUT_PDF = os.path.join(OUT_DIR, "Train_size_distribution.pdf")

MIN_LEN = 40
MAX_LEN = 500

MPE_COLOR = "#d62728"
BPE_COLOR = "#0D75BB"


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    df_len = pd.read_csv(LENGTH_MATRIX)
    df_group = pd.read_csv(GROUP_FILE)

    df_len["Sample"] = df_len["Sample"].astype(str).str.strip()
    df_group["Sample"] = df_group["Sample"].astype(str).str.strip()
    df_group["split"] = df_group["split"].astype(str).str.strip().str.lower()
    df_group["group"] = df_group["group"].astype(str).str.strip().str.lower()

    train_group = df_group[
        (df_group["split"] == "train")
        & (df_group["group"].isin(["mpe", "bpe"]))
    ][["Sample", "group"]].copy()

    train_group["display_group"] = train_group["group"].map({
        "mpe": "MPE",
        "bpe": "BPE",
    })

    df = df_len.merge(train_group, on="Sample", how="inner")

    length_pairs = []

    for col in df_len.columns:
        if col == "Sample":
            continue

        try:
            length = int(float(col))
        except ValueError:
            continue

        if MIN_LEN <= length <= MAX_LEN:
            length_pairs.append((col, length))

    length_pairs = sorted(length_pairs, key=lambda x: x[1])
    length_cols = [col for col, _ in length_pairs]
    length_values = np.array([length for _, length in length_pairs])

    df[length_cols] = df[length_cols].apply(pd.to_numeric, errors="coerce")

    summary_rows = []

    for display_group in ["BPE", "MPE"]:
        sub = df[df["display_group"] == display_group]

        mean_vals = sub[length_cols].mean(axis=0, skipna=True).to_numpy(dtype=float)
        sem_vals = sub[length_cols].sem(axis=0, skipna=True).to_numpy(dtype=float)

        for length, mean, sem in zip(length_values, mean_vals, sem_vals):
            summary_rows.append({
                "group": display_group,
                "length": length,
                "mean_relative_abundance": mean,
                "sem": sem,
            })

    summary_df = pd.DataFrame(summary_rows)

    plt.rcParams.update({
        "font.family": "Arial",
        "axes.linewidth": 1.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(8.5, 5.2), facecolor="white")
    ax.set_facecolor("white")

    for display_group, color in [("BPE", BPE_COLOR), ("MPE", MPE_COLOR)]:
        sub = summary_df[summary_df["group"] == display_group]

        x = sub["length"].to_numpy()
        y = sub["mean_relative_abundance"].to_numpy(dtype=float)
        sem = sub["sem"].to_numpy(dtype=float)

        ax.plot(
            x,
            y,
            color=color,
            linewidth=2.0,
            label=display_group,
            zorder=3,
        )

        ax.fill_between(
            x,
            y - sem,
            y + sem,
            color=color,
            alpha=0.18,
            linewidth=0,
            zorder=2,
        )

    ax.set_xlim(MIN_LEN, MAX_LEN)
    ax.set_xlabel("Fragment length (bp)", fontsize=13, fontweight="bold")
    ax.set_ylabel("Relative abundance", fontsize=13, fontweight="bold")

    ax.tick_params(axis="both", labelsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(
        frameon=False,
        fontsize=11,
        loc="upper right",
    )

    plt.tight_layout()
    fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight", format="pdf")
    plt.close(fig)



if __name__ == "__main__":
    main()
