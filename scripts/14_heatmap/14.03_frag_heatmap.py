"""Purpose: Frag heatmap for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial"]
mpl.rcParams["font.weight"] = "bold"
mpl.rcParams["axes.labelweight"] = "bold"
mpl.rcParams["axes.titleweight"] = "bold"

FEATURE_FILE = "results/14_heatmap/FRAG_feature.csv"
GROUP_FILE = "metadata/discovery_groupfile.csv"
OUT_DIR = "results/14_heatmap"
OUT_PDF = os.path.join(OUT_DIR, "Train_FRAG_subgroup_heatmap.pdf")

GROUP_ORDER = ["mpe", "bpe"]
SUBGROUP_ORDER = ["luad", "meso", "sclc", "scc", "tb", "ppe"]

group_label_map = {
    "mpe": "MPE",
    "bpe": "BPE",
}

group_colors = {
    "MPE": "#d62728",
    "BPE": "#0D75BB",
}

subgroup_colors = {
    "luad": "#D97C7C",
    "meso": "#B79AD6",
    "sclc": "#E6B65C",
    "scc": "#C07BA4",
    "tb": "#79B7D8",
    "ppe": "#7DBE9F",
}


def main():
    df_feat = pd.read_csv(FEATURE_FILE)
    df_group = pd.read_csv(GROUP_FILE)

    df_feat["Sample"] = df_feat["Sample"].astype(str).str.strip()

    df_group["Sample"] = df_group["Sample"].astype(str).str.strip()
    df_group["split"] = df_group["split"].astype(str).str.strip().str.lower()
    df_group["group"] = df_group["group"].astype(str).str.strip().str.lower()
    df_group["subgroup"] = df_group["subgroup"].astype(str).str.strip().str.lower()

    df_train = df_group[
        (df_group["split"] == "train") &
        (df_group["group"].isin(GROUP_ORDER)) &
        (df_group["subgroup"].isin(SUBGROUP_ORDER))
    ].copy()

    df_train["Group"] = df_train["group"].map(group_label_map)

    df_train["group_rank"] = pd.Categorical(
        df_train["group"],
        categories=GROUP_ORDER,
        ordered=True
    )

    df_train["subgroup_rank"] = pd.Categorical(
        df_train["subgroup"],
        categories=SUBGROUP_ORDER,
        ordered=True
    )

    df_train = df_train.sort_values(
        ["group_rank", "subgroup_rank", "Sample"]
    ).reset_index(drop=True)

    sample_cols = [
        s for s in df_train["Sample"].tolist()
        if s in df_feat["Sample"].values
    ]

    sample_info = df_train[df_train["Sample"].isin(sample_cols)].copy()
    sample_info = sample_info.set_index("Sample").loc[sample_cols].reset_index()

    df_feat = df_feat.set_index("Sample")
    feature_cols = df_feat.columns.tolist()

    mat = df_feat.loc[sample_cols, feature_cols].copy()
    mat = mat.apply(pd.to_numeric, errors="coerce")
    mat = mat.T

    mat_filled = mat.copy()
    row_median = mat_filled.median(axis=1)
    mat_filled = mat_filled.T.fillna(row_median).T

    col_median = mat_filled.median(axis=0)
    mat_filled = mat_filled.fillna(col_median)
    mat_filled = mat_filled.fillna(0)

    row_mean = mat_filled.mean(axis=1)
    row_std = mat_filled.std(axis=1, ddof=0).replace(0, pd.NA)

    mat_plot = mat_filled.sub(row_mean, axis=0).div(row_std, axis=0)
    mat_plot = mat_plot.fillna(0)

    row_linkage = linkage(
        pdist(mat_plot.values, metric="euclidean"),
        method="average"
    )

    cmap = LinearSegmentedColormap.from_list(
        "paper_blue_white_red",
        ["#86C5DA", "#F7F7F7", "#E64B35"],
        N=256
    )

    col_colors = pd.DataFrame({
        "Group": sample_info["Group"].map(group_colors).values,
        "Subgroup": sample_info["subgroup"].map(subgroup_colors).values,
    }, index=sample_info["Sample"])

    sns.set_theme(style="white", font_scale=0.75)

    g = sns.clustermap(
        mat_plot,
        row_cluster=True,
        col_cluster=False,
        row_linkage=row_linkage,
        col_colors=col_colors,
        cmap=cmap,
        vmin=-2,
        vmax=2,
        xticklabels=False,
        yticklabels=False,
        linewidths=0,
        figsize=(12, 10),
        dendrogram_ratio=(0.001, 0.001),
        cbar_kws={
            "label": "Z-score",
            "ticks": [-2, 0, 2],
        },
        cbar_pos=(0.98, 0.42, 0.018, 0.14),
        colors_ratio=(0.01, 0.014),
    )

    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)

    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("Significant fragment features")
    g.ax_heatmap.tick_params(axis="y", labelsize=7)

    g.ax_col_colors.set_ylabel("")
    g.ax_col_colors.axhline(
        1,
        color="white",
        linewidth=1.2,
        zorder=20
    )

    for spine in g.ax_col_colors.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.4)
        spine.set_color("#333333")

    g.ax_cbar.tick_params(labelsize=6, width=0.4, length=2)
    g.ax_cbar.set_ylabel(
        "Z-score",
        fontsize=9,
        fontweight="bold",
        labelpad=4
    )

    for spine in g.ax_cbar.spines.values():
        spine.set_linewidth(0.4)

    group_handles = [
        Patch(facecolor=group_colors["MPE"], edgecolor="none", label="MPE"),
        Patch(facecolor=group_colors["BPE"], edgecolor="none", label="BPE"),
    ]

    subgroup_handles = [
        Patch(facecolor=subgroup_colors[x], edgecolor="none", label=x.upper())
        for x in SUBGROUP_ORDER
    ]

    leg1 = g.ax_heatmap.legend(
        handles=group_handles,
        title="Group",
        loc="upper left",
        bbox_to_anchor=(1.04, 1.00),
        frameon=False,
        borderaxespad=0,
        fontsize=8,
        title_fontsize=8,
    )

    g.ax_heatmap.add_artist(leg1)

    g.ax_heatmap.legend(
        handles=subgroup_handles,
        title="Subgroup",
        loc="upper left",
        bbox_to_anchor=(1.04, 0.95),
        frameon=False,
        borderaxespad=0,
        fontsize=8,
        title_fontsize=8,
    )

    g.ax_heatmap.collections[0].set_rasterized(True)

    plt.savefig(
        OUT_PDF,
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
        transparent=False
    )
    plt.show()


if __name__ == "__main__":
    main()
