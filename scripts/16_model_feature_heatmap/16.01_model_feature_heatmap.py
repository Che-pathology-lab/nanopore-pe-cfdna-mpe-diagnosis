"""Purpose: Model feature heatmap for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist
import matplotlib as mpl

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial"]
mpl.rcParams["font.weight"] = "bold"
mpl.rcParams["axes.labelweight"] = "bold"
mpl.rcParams["axes.titleweight"] = "bold"

feature_file = "results/05_nanopore_multimodal/discovery_Multi_features_wide_table.csv"
group_file = "metadata/discovery_groupfile.csv"

meth_file = "results/02_nanopore_methylation/discovery_METH_prediction_results_with_scores.csv"
cna_file = "results/03_nanopore_cna/discovery_CNA_prediction_results_with_scores.csv"
frag_file = "results/04_nanopore_fragmentomics/discovery_FRAG_prediction_results_with_scores.csv"
multi_file = "results/05_nanopore_multimodal/discovery_Multi_prediction_results_with_scores.csv"
output_pdf = "results/16_model_feature_heatmap/train_features_and_model_labels_heatmap.pdf"
output_order_csv = "results/16_model_feature_heatmap/train_features_heatmap_order.csv"


features = [
    "Global_Hyper_Mean",
    "Global_Hypo_Mean",
    "Max_Oscillating_Length",
    "Mean_Changepoint_Jump",
    "Mean_Aberrant_Seg_Size",
    "AACT",
    "CGAT",
    "CTCG",
    "CTGT",
    "CTTG",
]

group_label_map = {
    "mpe": "MPE",
    "bpe": "BPE",
}

group_colors = {
    "MPE": "#d62728",
    "BPE": "#0D75BB",
}


def read_pred(file, label_col, out_col, prob_col, out_prob_col):
    df = pd.read_csv(file)
    df["Sample"] = df["Sample"].astype(str).str.strip()

    df[out_col] = pd.to_numeric(df[label_col], errors="coerce")
    df[out_prob_col] = pd.to_numeric(df[prob_col], errors="coerce")

    return df[["Sample", out_col, out_prob_col]]


def main() -> None:
    df_feat = pd.read_csv(feature_file)
    df_group = pd.read_csv(group_file)

    df_feat["Sample"] = df_feat["Sample"].astype(str).str.strip()
    df_group["Sample"] = df_group["Sample"].astype(str).str.strip()

    df_group_train = df_group[
        (df_group["split"].astype(str).str.lower() == "train")
        & (df_group["group"].astype(str).str.lower().isin(["mpe", "bpe"]))
    ].copy()

    df_group_train["group_lower"] = df_group_train["group"].astype(str).str.lower()
    df_group_train["Group"] = df_group_train["group_lower"].map(group_label_map)

    df_feature = df_group_train[["Sample", "Group"]].merge(
        df_feat[["Sample"] + features],
        on="Sample",
        how="inner",
    )

    df_meth = read_pred(
        meth_file,
        label_col="METH_Predicted_Label",
        out_col="METH",
        prob_col="METH_Score_Prob",
        out_prob_col="METH_Score_Prob",
    )

    df_cna = read_pred(
        cna_file,
        label_col="CNA_Predicted_Label",
        out_col="CNA",
        prob_col="CNA_Score_Prob",
        out_prob_col="CNA_Score_Prob",
    )

    df_frag = read_pred(
        frag_file,
        label_col="FRAG_Predicted_Label",
        out_col="FRAG",
        prob_col="FRAG_Score_Prob",
        out_prob_col="FRAG_Score_Prob",
    )

    df_multi = read_pred(
        multi_file,
        label_col="Multi_Predicted_Label",
        out_col="Multi",
        prob_col="Multi_Score_Prob",
        out_prob_col="Multi_Score_Prob",
    )

    df_label = df_group_train[["Sample", "Group"]].copy()
    df_label = df_label.merge(df_meth, on="Sample", how="inner")
    df_label = df_label.merge(df_cna, on="Sample", how="inner")
    df_label = df_label.merge(df_frag, on="Sample", how="inner")
    df_label = df_label.merge(df_multi, on="Sample", how="inner")

    common_samples = sorted(set(df_feature["Sample"]) & set(df_label["Sample"]))

    df_feature = df_feature[df_feature["Sample"].isin(common_samples)].copy()
    df_label = df_label[df_label["Sample"].isin(common_samples)].copy()

    df_feature = df_feature.sort_values(["Group", "Sample"], ascending=[False, True])

    mat = df_feature.set_index("Sample")[features].apply(pd.to_numeric, errors="coerce")

    std = mat.std(axis=0, ddof=0).replace(0, np.nan)
    mat_z = (mat - mat.mean(axis=0)) / std
    mat_z = mat_z.fillna(0)

    plot_feature = mat_z.T

    col_linkage = linkage(
        pdist(plot_feature.T, metric="euclidean"),
        method="average",
    )

    ordered_samples = plot_feature.columns[leaves_list(col_linkage)].tolist()
    plot_feature = plot_feature[ordered_samples]

    model_cols = [
        "METH",
        "CNA",
        "FRAG",
        "Multi",
    ]

    prob_cols = [
        "METH_Score_Prob",
        "CNA_Score_Prob",
        "FRAG_Score_Prob",
        "Multi_Score_Prob",
    ]

    model_label_map = {
        "METH": "METH",
        "CNA": "CNA",
        "FRAG": "FRAG",
        "Multi": "Multi",
    }

    order_df = pd.DataFrame({
        "heatmap_order": np.arange(1, len(ordered_samples) + 1),
        "Sample": ordered_samples,
    })

    order_df = order_df.merge(
        df_feature[["Sample", "Group"]],
        on="Sample",
        how="left",
    )

    order_df = order_df.merge(
        df_label[["Sample"] + model_cols + prob_cols],
        on="Sample",
        how="left",
    )

    order_df["True_Label_Numeric"] = order_df["Group"].map({
        "BPE": 0,
        "MPE": 1,
    })

    order_df["METH_Label"] = order_df["METH"].map({0: "BPE", 1: "MPE"})
    order_df["CNA_Label"] = order_df["CNA"].map({0: "BPE", 1: "MPE"})
    order_df["FRAG_Label"] = order_df["FRAG"].map({0: "BPE", 1: "MPE"})
    order_df["Multi_Label"] = order_df["Multi"].map({0: "BPE", 1: "MPE"})

    order_df.to_csv(output_order_csv, index=False)
    df_label = df_label.set_index("Sample").loc[ordered_samples].reset_index()
    plot_label = df_label.set_index("Sample")[model_cols].T
    plot_label.index = [model_label_map[x] for x in plot_label.index]
    sample_group = df_feature.set_index("Sample").loc[ordered_samples, "Group"]

    group_numeric = pd.DataFrame(
        [sample_group.map({"BPE": 0, "MPE": 1}).values],
        columns=ordered_samples,
        index=["True label"],
    )

    feature_cmap = LinearSegmentedColormap.from_list(
        "paper_blue_white_red",
        ["#86C5DA", "#F7F7F7", "#E64B35"],
        N=256,
    )

    label_cmap = ListedColormap(["#86C5DA", "#E64B35"])
    label_norm = BoundaryNorm([-0.5, 0.5, 1.5], label_cmap.N)

    group_cmap = ListedColormap([group_colors["BPE"], group_colors["MPE"]])
    group_norm = BoundaryNorm([-0.5, 0.5, 1.5], group_cmap.N)

    sns.set_theme(style="white", font_scale=0.8)

    fig = plt.figure(figsize=(13, 7.2))

    gs = fig.add_gridspec(
        nrows=4,
        ncols=2,
        width_ratios=[40, 0.7],
        height_ratios=[0.55, 0.18, 2.2, 0.8],
        wspace=0.08,
        hspace=0.04,
    )

    ax_dend = fig.add_subplot(gs[0, 0])
    ax_group = fig.add_subplot(gs[1, 0])
    ax_feature = fig.add_subplot(gs[2, 0])
    ax_label = fig.add_subplot(gs[3, 0])

    cax_feature = fig.add_subplot(gs[2, 1])
    cax_label = fig.add_subplot(gs[3, 1])

    dendrogram(
        col_linkage,
        ax=ax_dend,
        no_labels=True,
        color_threshold=0,
        above_threshold_color="#666666",
    )

    for line in ax_dend.collections:
        line.set_linewidth(0.45)
        line.set_color("#666666")

    ax_dend.set_xticks([])
    ax_dend.set_yticks([])

    for spine in ax_dend.spines.values():
        spine.set_visible(False)

    sns.heatmap(
        group_numeric,
        ax=ax_group,
        cmap=group_cmap,
        norm=group_norm,
        cbar=False,
        xticklabels=False,
        yticklabels=True,
        linewidths=0,
    )

    ax_group.set_xlabel("")
    ax_group.set_ylabel("")
    ax_group.tick_params(axis="y", labelrotation=0, labelsize=8, length=0)

    sns.heatmap(
        plot_feature,
        ax=ax_feature,
        cmap=feature_cmap,
        vmin=-2,
        vmax=2,
        center=0,
        cbar=True,
        linewidths=0.25,
        linecolor="#D0D0D0",
        cbar_ax=cax_feature,
        cbar_kws={"ticks": [-2, 0, 2]},
        xticklabels=False,
        yticklabels=True,
    )

    ax_feature.set_xlabel("")
    ax_feature.set_ylabel("")
    ax_feature.tick_params(axis="y", labelsize=8)

    cax_feature.set_ylabel("Z-score", fontsize=8, labelpad=2)
    cax_feature.tick_params(labelsize=8, width=0.3, length=1.8)

    sns.heatmap(
        plot_label,
        ax=ax_label,
        cmap=label_cmap,
        norm=label_norm,
        cbar=True,
        cbar_ax=cax_label,
        cbar_kws={"ticks": [0, 1]},
        xticklabels=False,
        yticklabels=True,
        linewidths=0.25,
        linecolor="#D0D0D0",
    )

    ax_label.set_xlabel("")
    ax_label.set_ylabel("")
    ax_label.tick_params(axis="y", labelsize=8)

    ax_label.set_yticklabels(
        ax_label.get_yticklabels(),
        rotation=0,
        ha="right",
        va="center",
        fontsize=8,
        fontweight="bold",
    )

    cax_label.set_yticklabels(["BPE", "MPE"], fontsize=8)
    cax_label.set_ylabel("Model Prediction", fontsize=8, labelpad=3)
    cax_label.tick_params(labelsize=8, width=0.3, length=1.8)

    for ax in [cax_feature, cax_label]:
        for spine in ax.spines.values():
            spine.set_linewidth(0.4)

    legend_handles = [
        Patch(facecolor=group_colors["MPE"], edgecolor="none", label="MPE"),
        Patch(facecolor=group_colors["BPE"], edgecolor="none", label="BPE"),
    ]

    leg = ax_feature.legend(
        handles=legend_handles,
        title="True label",
        loc="upper left",
        bbox_to_anchor=(1.08, 1.03),
        frameon=False,
        borderaxespad=0,
        fontsize=8,
        title_fontsize=8,
    )

    leg.get_title().set_fontweight("bold")

    for text in leg.get_texts():
        text.set_fontweight("bold")

    plt.savefig(output_pdf, bbox_inches="tight")
    plt.show()



if __name__ == "__main__":
    main()
