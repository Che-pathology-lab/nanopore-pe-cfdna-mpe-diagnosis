"""Purpose: Multi downsampling model for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import joblib
import re
from sklearn.metrics import roc_auc_score, confusion_matrix


MODEL_FILE = "results/05_nanopore_multimodal/discovery_Multi_logistic_model.joblib"
GROUP_FILE = "metadata/discovery_groupfile.csv"
DEPTHS = ("0.01", "0.1")
FEATURE_FILE_TEMPLATE = (
    "results/06_downsampling/Multi/"
    "discovery_{depth}_Multi_features_wide_table.csv"
)
PREDICTION_CSV_TEMPLATE = (
    "results/06_downsampling/Multi/"
    "downsampling_{depth}x_Multi_prediction_results_with_scores.csv"
)

def run_depth(depth):

    model_pkg = joblib.load(MODEL_FILE)
    pipeline = model_pkg["pipeline"]
    threshold = model_pkg["threshold"]
    required_features = model_pkg["features"]

    feature_file = FEATURE_FILE_TEMPLATE.format(depth=depth)
    prediction_csv = PREDICTION_CSV_TEMPLATE.format(depth=depth)

    df_features = pd.read_csv(feature_file)
    df_group = pd.read_csv(GROUP_FILE)
    df_group["Sample"] = df_group["Sample"].astype(str).str.strip()
    df_group["group"] = df_group["group"].astype(str).str.strip().str.lower()
    df_group["split"] = df_group["split"].astype(str).str.strip().str.lower()

    df_features["Original_Sample"] = df_features["Sample"].apply(
        lambda x: re.sub(r"_\d+\.?\d*x_rep\d+$", "", str(x))
    )

    df_features["Replicate"] = df_features["Sample"].apply(
        lambda x: re.search(r"(rep\d+)", str(x)).group(1)
        if re.search(r"(rep\d+)", str(x))
        else "unknown"
    )

    validation_samples = set(
        df_group.loc[
            df_group["split"] == "test",
            "Sample",
        ]
    )
    if not validation_samples:
        raise ValueError(
            f"No samples with split='test' were found in {GROUP_FILE}"
        )

    df_features = df_features[
        df_features["Original_Sample"].isin(validation_samples)
    ].copy()

    group_map = dict(zip(df_group["Sample"], df_group["group"]))
    df_features["group"] = df_features["Original_Sample"].map(group_map)

    label_mapping = {"mpe": 1, "bpe": 0}

    df_full = df_features.copy()
    df_full = df_full[df_full["group"].isin(label_mapping.keys())].copy()
    df_full.sort_values(by="Sample", inplace=True)
    df_full.reset_index(drop=True, inplace=True)

    observed_samples = set(df_full["Original_Sample"])
    if observed_samples != validation_samples:
        missing = sorted(validation_samples - observed_samples)
        extra = sorted(observed_samples - validation_samples)
        raise ValueError(
            f"Validation sample mismatch. Missing={missing}, extra={extra}"
        )

    if df_full["Sample"].duplicated().any():
        duplicates = sorted(
            df_full.loc[df_full["Sample"].duplicated(), "Sample"].unique()
        )
        raise ValueError(
            f"Duplicate downsampling sample rows found: {duplicates}"
        )

    missing_features = [f for f in required_features if f not in df_full.columns]
    if missing_features:
        raise ValueError(f"Missing required features: {missing_features}")

    X_all = df_full[required_features]
    y_all = df_full["group"].map(label_mapping)

    y_pred_prob = pipeline.predict_proba(X_all)[:, 1]
    y_pred_logit = pipeline.decision_function(X_all)
    y_pred_label = (y_pred_prob >= threshold).astype(int)

    cm = confusion_matrix(y_all, y_pred_label, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()

    auc_test = roc_auc_score(y_all, y_pred_prob)
    print(f"Downsampling Multi {depth}x AUC: {auc_test:.4f}")
    print(f"Confusion Matrix: TN={tn}, FP={fp}, FN={fn}, TP={tp}")

    res_df = df_full[["Sample", "Original_Sample", "Replicate"]].copy()
    res_df["Multi_Score_Prob"] = y_pred_prob
    res_df["Multi_Score_Logit"] = y_pred_logit
    res_df["Multi_Predicted_Label"] = y_pred_label

    columns_to_save = [
        "Sample",
        "Original_Sample",
        "Replicate",
        "Multi_Score_Prob",
        "Multi_Score_Logit",
        "Multi_Predicted_Label",
    ]

    res_df[columns_to_save].to_csv(prediction_csv, index=False)


def main():
    for depth in DEPTHS:
        run_depth(depth)


if __name__ == "__main__":
    main()
