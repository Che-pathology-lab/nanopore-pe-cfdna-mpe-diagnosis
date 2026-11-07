"""Purpose: Measure tumor-fraction agreement across completed PoN analyses.

Input: A table with Sample and at least two columns named ``<pon_name>_TF``,
       plus discovery metadata.
Outputs: Pearson matrix, pairwise/summary tables, a boxplot, and a heatmap.
Usage: python scripts/18_pon_robustness/18.02_tumor_fraction_correlation.py
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


TUMOR_FRACTION_FILE = "results/18_pon_robustness/pon_tumor_fraction.csv"
GROUP_FILE = "metadata/discovery_groupfile.csv"
OUTPUT_DIR = "results/18_pon_robustness"

CORRELATION_OUTPUT = os.path.join(OUTPUT_DIR, "tumor_fraction_correlation_matrix.csv")
PAIRWISE_OUTPUT = os.path.join(OUTPUT_DIR, "tumor_fraction_pairwise_correlations.csv")
SUMMARY_OUTPUT = os.path.join(OUTPUT_DIR, "tumor_fraction_correlation_summary.csv")
BOXPLOT_OUTPUT = os.path.join(OUTPUT_DIR, "tumor_fraction_correlation_boxplot.pdf")
HEATMAP_OUTPUT = os.path.join(OUTPUT_DIR, "tumor_fraction_correlation_heatmap.pdf")

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "Arial"


def main() -> None:
    tumor_fraction = pd.read_csv(TUMOR_FRACTION_FILE)
    metadata = pd.read_csv(GROUP_FILE)
    tumor_fraction["Sample"] = tumor_fraction["Sample"].astype(str).str.strip()
    metadata["Sample"] = metadata["Sample"].astype(str).str.strip()
    metadata["split"] = metadata["split"].astype(str).str.strip().str.lower()

    train_samples = metadata.loc[metadata["split"] == "train", ["Sample"]]
    data = train_samples.merge(tumor_fraction, on="Sample", how="inner")
    tumor_fraction_columns = [
        column for column in data.columns if column.endswith("_TF")
    ]
    data[tumor_fraction_columns] = data[tumor_fraction_columns].apply(
        pd.to_numeric, errors="coerce"
    )

    correlation = data[tumor_fraction_columns].corr(method="pearson")
    correlation.to_csv(CORRELATION_OUTPUT)

    upper_triangle = np.triu(np.ones(correlation.shape), k=1).astype(bool)
    pairwise = correlation.where(upper_triangle).stack().reset_index()
    pairwise.columns = ["PoN_1", "PoN_2", "pearson_r"]
    pairwise.to_csv(PAIRWISE_OUTPUT, index=False)

    summary = pd.DataFrame(
        [
            {
                "n_TF_columns": len(tumor_fraction_columns),
                "n_pairwise_comparisons": len(pairwise),
                "mean_pearson_r": pairwise["pearson_r"].mean(),
                "median_pearson_r": pairwise["pearson_r"].median(),
                "sd_pearson_r": pairwise["pearson_r"].std(),
                "q25_pearson_r": pairwise["pearson_r"].quantile(0.25),
                "q75_pearson_r": pairwise["pearson_r"].quantile(0.75),
                "min_pearson_r": pairwise["pearson_r"].min(),
                "max_pearson_r": pairwise["pearson_r"].max(),
            }
        ]
    )
    summary.to_csv(SUMMARY_OUTPUT, index=False)

    figure, axis = plt.subplots(figsize=(4.2, 5.2))
    sns.boxplot(
        y=pairwise["pearson_r"],
        color="#9ecae1",
        width=0.35,
        showfliers=False,
        ax=axis,
    )
    sns.stripplot(
        y=pairwise["pearson_r"],
        color="black",
        alpha=0.45,
        size=4,
        jitter=0.15,
        ax=axis,
    )
    axis.set_ylim(0, 1.12)
    axis.set_xlabel("Tumor fraction")
    axis.set_ylabel("Pairwise Pearson correlation")
    sns.despine()
    figure.tight_layout()
    figure.savefig(BOXPLOT_OUTPUT, dpi=300, bbox_inches="tight")
    plt.close(figure)

    heatmap_size = max(6, 0.45 * len(tumor_fraction_columns))
    figure, axis = plt.subplots(figsize=(heatmap_size, heatmap_size))
    sns.heatmap(
        correlation,
        cmap="vlag",
        vmin=-1,
        vmax=1,
        center=0,
        square=True,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Pearson correlation", "shrink": 0.75},
        ax=axis,
    )
    axis.set_title("Tumor-fraction correlation")
    figure.tight_layout()
    figure.savefig(HEATMAP_OUTPUT, dpi=300, bbox_inches="tight")
    plt.close(figure)

    print(f"Saved: {CORRELATION_OUTPUT}")
    print(f"Saved: {PAIRWISE_OUTPUT}")
    print(f"Saved: {SUMMARY_OUTPUT}")
    print(f"Saved: {BOXPLOT_OUTPUT}")
    print(f"Saved: {HEATMAP_OUTPUT}")


if __name__ == "__main__":
    main()
