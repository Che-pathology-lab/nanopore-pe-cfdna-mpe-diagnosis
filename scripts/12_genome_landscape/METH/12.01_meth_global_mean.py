import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import mannwhitneyu

METH_FILE = "results/02_nanopore_methylation/threshold_10_discovery_cohort_methylationfile.tsv.gz"
GROUP_FILE = "metadata/discovery_groupfile.csv"
OUT_PDF = "results/12_genome_landscape/METH/Train_METH_genome_wide_mean.pdf"

PALETTE = {"BPE": "#0D75BB", "MPE": "#d62728"}
GROUP_MAP = {"bpe": "BPE", "mpe": "MPE"}
GROUP_ORDER = ["BPE", "MPE"]


def get_significance_label(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def get_axis_limits(values):
    values = values.dropna().astype(float)
    data_min = values.min()
    data_max = values.max()
    data_range = max(data_max - data_min, 2.0)
    lower = np.floor(data_min - 0.12 * data_range)
    upper = np.ceil(data_max + 0.35 * data_range)
    span = upper - lower
    step = 1 if span <= 8 else 2 if span <= 16 else 5
    ticks = np.arange(np.floor(lower / step) * step, np.ceil(upper / step) * step + step, step)
    return (lower, upper), ticks


def main():
    plt.rcParams.update({
        "font.sans-serif": ["Arial"],
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 2.0,
    })
    sns.set_theme(style="ticks")

    df = pd.read_csv(METH_FILE, sep="\t")
    df.columns = [str(column).strip() for column in df.columns]
    df = df[df.iloc[:, 0].isin([f"chr{i}" for i in range(1, 23)])].copy()

    meth_df = df.drop(columns=df.columns[:3]).apply(pd.to_numeric, errors="coerce")

    group_df = pd.read_csv(GROUP_FILE)
    group_df["Sample"] = group_df["Sample"].astype(str).str.strip()
    group_df["split"] = group_df["split"].astype(str).str.strip().str.lower()
    group_df["group"] = group_df["group"].astype(str).str.strip().str.lower()

    train_df = group_df[
        (group_df["split"] == "train")
        & group_df["group"].isin(["mpe", "bpe"])
        & group_df["Sample"].isin(meth_df.columns)
    ].copy()

    train_samples = train_df["Sample"].tolist()
    meth_mat = meth_df[train_samples]
    callable_mask = meth_mat.notna().mean(axis=1) >= 0.90
    sample_mean = meth_mat.loc[callable_mask].mean(axis=0, skipna=True)

    summary_df = train_df[["Sample", "group"]].copy()
    summary_df["Mean_Methylation"] = summary_df["Sample"].map(sample_mean)
    summary_df = summary_df.dropna(subset=["Mean_Methylation"])

    if summary_df["Mean_Methylation"].max() <= 1.5:
        summary_df["Mean_Methylation"] *= 100

    summary_df["group_label"] = summary_df["group"].map(GROUP_MAP)

    bpe = summary_df.loc[summary_df["group"] == "bpe", "Mean_Methylation"]
    mpe = summary_df.loc[summary_df["group"] == "mpe", "Mean_Methylation"]
    _, p_value = mannwhitneyu(mpe, bpe, alternative="two-sided")
    significance = get_significance_label(p_value)

    fig, ax = plt.subplots(figsize=(5, 6))
    sns.violinplot(
        data=summary_df,
        x="group_label",
        y="Mean_Methylation",
        hue="group_label",
        order=GROUP_ORDER,
        hue_order=GROUP_ORDER,
        palette=PALETTE,
        linewidth=1.5,
        inner=None,
        cut=2,
        saturation=0.9,
        legend=False,
        ax=ax,
    )

    box_data = [
        summary_df.loc[summary_df["group_label"] == group, "Mean_Methylation"].dropna()
        for group in GROUP_ORDER
    ]
    ax.boxplot(
        box_data,
        positions=[0, 1],
        widths=0.12,
        patch_artist=True,
        showfliers=False,
        manage_ticks=False,
        boxprops=dict(facecolor="white", edgecolor="black", linewidth=1, alpha=0.75),
        medianprops=dict(color="black", linewidth=1),
        whiskerprops=dict(color="black", linewidth=1),
        capprops=dict(color="black", linewidth=1),
    )

    ylim, yticks = get_axis_limits(summary_df["Mean_Methylation"])
    ax.set_ylim(ylim)
    ax.set_yticks(yticks)

    if significance != "ns":
        span = ylim[1] - ylim[0]
        y = ylim[1] - 0.08 * span
        ax.plot([0, 0, 1, 1], [y - 0.05 * span, y, y, y - 0.05 * span], color="black", lw=1.5)
        ax.text(0.5, y + 0.02 * span, significance, ha="center", va="bottom", fontsize=15, fontweight="bold")

    ax.set_xlabel("")
    ax.set_ylabel("Genome-wide mean methylation (%)", fontsize=13, fontweight="bold")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(GROUP_ORDER, fontsize=13, fontweight="bold")
    ax.tick_params(axis="x", width=1.5, length=5)
    ax.tick_params(axis="y", labelsize=12, width=1.5, length=5)
    ax.legend(
        handles=[mpatches.Patch(facecolor=PALETTE[group], label=group) for group in GROUP_ORDER],
        loc="upper right",
        frameon=False,
        fontsize=11,
    )

    sns.despine()
    fig.subplots_adjust(left=0.22, right=0.95, bottom=0.15, top=0.92)
    fig.savefig(OUT_PDF, dpi=300)


if __name__ == "__main__":
    main()
