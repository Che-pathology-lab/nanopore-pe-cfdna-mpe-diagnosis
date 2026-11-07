import pandas as pd
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "Arial"

MODEL_FILE = "results/09_ulp_validation/ULP_Multi_prediction_results_with_scores.csv"
GROUP_FILE = "metadata/ULP_groupfile.csv"
OUTPUT_DIR = "results/10_cytology_comparison"
OUTPUT_PDF = os.path.join(OUTPUT_DIR, "ULP_cytology_subtype_comparison.pdf")

SUBGROUPS = ["LUAD", "SCLC", "MESO", "SCC"]
METHODS = ["Combined Cytology", "Multimodal Model", "CURE-guided Cytology"]
COLORS = {
    "Combined Cytology": "#0D75BB",
    "Multimodal Model": "#D62728",
    "CURE-guided Cytology": "#E69F00",
}


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df_model = pd.read_csv(MODEL_FILE)
    df_group = pd.read_csv(GROUP_FILE)

    for df in [df_model, df_group]:
        df["Sample"] = df["Sample"].astype(str).str.strip()

    df_group["split"] = df_group["split"].astype(str).str.strip().str.lower()
    df_group["group"] = df_group["group"].astype(str).str.strip().str.lower()
    df_group["subgroup"] = df_group["subgroup"].astype(str).str.strip().str.upper()

    df_merge = pd.merge(
        df_group[["Sample", "group", "subgroup", "split", "Combined_Cyto"]],
        df_model[["Sample", "Multi_Predicted_Label"]],
        on="Sample",
        how="inner",
    )

    df_merge = df_merge[
        (df_merge["split"] == "ulp") & df_merge["group"].isin(["mpe", "bpe"])
    ].copy()
    df_merge["y_true"] = df_merge["group"].map({"mpe": 1, "bpe": 0})
    df_merge["Combined_Cyto"] = pd.to_numeric(df_merge["Combined_Cyto"], errors="coerce")
    df_merge["Multi_Predicted_Label"] = pd.to_numeric(
        df_merge["Multi_Predicted_Label"], errors="coerce"
    )
    df_merge = df_merge.dropna(subset=["Combined_Cyto", "Multi_Predicted_Label"])
    df_merge["CURE_Guided_Pred"] = (
        (df_merge["Combined_Cyto"] == 1) | (df_merge["Multi_Predicted_Label"] == 1)
    ).astype(int)

    subgroup_results = []
    for subgroup in SUBGROUPS:
        df_sub = df_merge[(df_merge["subgroup"] == subgroup) & (df_merge["y_true"] == 1)]
        n_mpe = len(df_sub)

        subgroup_results.append({
            "Subgroup": subgroup,
            "N_Samples": n_mpe,
            "Combined Cytology": (df_sub["Combined_Cyto"] == 1).sum() / n_mpe * 100,
            "Multimodal Model": (df_sub["Multi_Predicted_Label"] == 1).sum() / n_mpe * 100,
            "CURE-guided Cytology": (df_sub["CURE_Guided_Pred"] == 1).sum() / n_mpe * 100,
        })

    subgroup_df = pd.DataFrame(subgroup_results)
    plot_data = subgroup_df.melt(
        id_vars=["Subgroup", "N_Samples"],
        value_vars=METHODS,
        var_name="Method",
        value_name="Sensitivity",
    )
    plot_data["Subgroup"] = plot_data.apply(
        lambda row: f'{row["Subgroup"]}\n(n={int(row["N_Samples"])})', axis=1
    )

    sns.set_theme(style="ticks", rc={"font.family": "Arial"})
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    sns.barplot(
        data=plot_data,
        x="Subgroup",
        y="Sensitivity",
        hue="Method",
        hue_order=METHODS,
        palette=COLORS,
        saturation=0.85,
        ax=ax,
    )

    for container in ax.containers:
        ax.bar_label(container, fmt="%.1f", padding=4, fontsize=8, fontweight="bold")

    ax.set_ylim(0, 115)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_xlabel("")
    ax.set_ylabel("Sensitivity (%)", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", labelsize=10)
    ax.tick_params(axis="y", labelsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.14), ncol=3, frameon=False, fontsize=8)
    sns.despine()

    fig.tight_layout()
    fig.savefig(OUTPUT_PDF, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
