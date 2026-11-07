"""Purpose: Measure CNA-profile agreement across completed PoN analyses.

Input: A wide CNA table with columns named
       ``<pon_name>__<Sample>__Corrected_Copy_Number`` and discovery metadata.
Outputs: Pairwise correlations, sample/overall summaries, and a boxplot.
Usage: python scripts/18_pon_robustness/18.01_cna_profile_correlation.py
"""

import itertools
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr


CNA_FILE = "results/18_pon_robustness/pon_Corrected_Copy_Number.csv"
GROUP_FILE = "metadata/discovery_groupfile.csv"
OUTPUT_DIR = "results/18_pon_robustness"

PAIRWISE_OUTPUT = os.path.join(OUTPUT_DIR, "CNA_profile_pairwise_correlations.csv")
SAMPLE_OUTPUT = os.path.join(OUTPUT_DIR, "CNA_profile_sample_summary.csv")
SUMMARY_OUTPUT = os.path.join(OUTPUT_DIR, "CNA_profile_overall_summary.csv")
PLOT_OUTPUT = os.path.join(OUTPUT_DIR, "CNA_profile_pairwise_pearson_boxplot.pdf")

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "Arial"


def parse_cna_column(column: str) -> tuple[str, str] | None:
    suffix = "__Corrected_Copy_Number"
    if not column.endswith(suffix):
        return None

    parts = column[: -len(suffix)].split("__")
    if len(parts) != 2:
        return None
    return parts[0], parts[1]


def safe_correlation(x: pd.Series, y: pd.Series, method: str) -> float:
    valid = x.notna() & y.notna()
    x_values = x[valid].astype(float).to_numpy()
    y_values = y[valid].astype(float).to_numpy()

    if len(x_values) < 3 or np.std(x_values) == 0 or np.std(y_values) == 0:
        return np.nan
    if method == "pearson":
        return pearsonr(x_values, y_values)[0]
    return spearmanr(x_values, y_values)[0]


def main() -> None:
    cna = pd.read_csv(CNA_FILE)
    metadata = pd.read_csv(GROUP_FILE)
    metadata["Sample"] = metadata["Sample"].astype(str).str.strip()
    metadata["split"] = metadata["split"].astype(str).str.strip().str.lower()
    train_samples = set(metadata.loc[metadata["split"] == "train", "Sample"])

    column_information = []
    for column in cna.columns:
        parsed = parse_cna_column(column)
        if parsed is not None and parsed[1] in train_samples:
            column_information.append(
                {"column": column, "pon_name": parsed[0], "Sample": parsed[1]}
            )

    information = pd.DataFrame(column_information)
    cna[information["column"]] = cna[information["column"]].apply(
        pd.to_numeric, errors="coerce"
    )

    pairwise_rows = []
    for sample, sample_information in information.groupby("Sample"):
        records = sample_information.to_dict("records")
        for first, second in itertools.combinations(records, 2):
            x = cna[first["column"]]
            y = cna[second["column"]]
            valid = x.notna() & y.notna()
            difference = x[valid].astype(float).to_numpy() - y[valid].astype(float).to_numpy()
            pairwise_rows.append(
                {
                    "Sample": sample,
                    "PoN_1": first["pon_name"],
                    "PoN_2": second["pon_name"],
                    "n_bins_compared": int(valid.sum()),
                    "pearson_r": safe_correlation(x, y, "pearson"),
                    "spearman_r": safe_correlation(x, y, "spearman"),
                    "mae": np.mean(np.abs(difference)) if len(difference) else np.nan,
                    "rmse": np.sqrt(np.mean(difference**2)) if len(difference) else np.nan,
                }
            )

    pairwise = pd.DataFrame(pairwise_rows)
    pairwise.to_csv(PAIRWISE_OUTPUT, index=False)

    sample_summary = (
        pairwise.groupby("Sample")
        .agg(
            n_pairwise=("pearson_r", "count"),
            median_pearson_r=("pearson_r", "median"),
            mean_pearson_r=("pearson_r", "mean"),
            min_pearson_r=("pearson_r", "min"),
            q25_pearson_r=("pearson_r", lambda values: values.quantile(0.25)),
            q75_pearson_r=("pearson_r", lambda values: values.quantile(0.75)),
            median_spearman_r=("spearman_r", "median"),
            median_mae=("mae", "median"),
            median_rmse=("rmse", "median"),
        )
        .reset_index()
    )
    sample_summary.to_csv(SAMPLE_OUTPUT, index=False)

    overall_summary = pd.DataFrame(
        [
            {
                "n_samples": sample_summary["Sample"].nunique(),
                "n_pairwise_comparisons": len(pairwise),
                "mean_pearson_r": pairwise["pearson_r"].mean(),
                "median_pearson_r": pairwise["pearson_r"].median(),
                "q25_pearson_r": pairwise["pearson_r"].quantile(0.25),
                "q75_pearson_r": pairwise["pearson_r"].quantile(0.75),
                "min_pearson_r": pairwise["pearson_r"].min(),
                "max_pearson_r": pairwise["pearson_r"].max(),
                "median_spearman_r": pairwise["spearman_r"].median(),
                "median_mae": pairwise["mae"].median(),
                "median_rmse": pairwise["rmse"].median(),
            }
        ]
    )
    overall_summary.to_csv(SUMMARY_OUTPUT, index=False)

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
        alpha=0.35,
        size=3,
        jitter=0.18,
        ax=axis,
    )
    axis.set_ylim(0, 1.12)
    axis.set_xlabel("CNA profile")
    axis.set_ylabel("Pairwise Pearson correlation")
    sns.despine()
    figure.tight_layout()
    figure.savefig(PLOT_OUTPUT, dpi=300, bbox_inches="tight")
    plt.close(figure)

    print(f"Saved: {PAIRWISE_OUTPUT}")
    print(f"Saved: {SAMPLE_OUTPUT}")
    print(f"Saved: {SUMMARY_OUTPUT}")
    print(f"Saved: {PLOT_OUTPUT}")


if __name__ == "__main__":
    main()
