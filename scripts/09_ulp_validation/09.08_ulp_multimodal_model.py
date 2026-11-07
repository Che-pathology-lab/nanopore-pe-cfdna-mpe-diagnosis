"""Purpose: Apply the discovery Multi model to the ULP cohort.

Edit the paths below, then run this script directly.
"""

import joblib
import pandas as pd
from sklearn.metrics import confusion_matrix, roc_auc_score


MODEL_FILE = "results/05_nanopore_multimodal/discovery_Multi_logistic_model.joblib"
FEATURE_FILE = "results/09_ulp_validation/ULP_Multi_features.csv"
GROUP_FILE = "metadata/ULP_groupfile.csv"
PREDICTION_CSV = "results/09_ulp_validation/ULP_Multi_prediction_results_with_scores.csv"


def main() -> None:
    model_package = joblib.load(MODEL_FILE)
    pipeline = model_package["pipeline"]
    threshold = model_package["threshold"]
    required_features = model_package["features"]

    df_features = pd.read_csv(FEATURE_FILE)
    group_df = pd.read_csv(GROUP_FILE)
    group_map = dict(zip(group_df["Sample"], group_df["group"]))
    df_features["group"] = df_features["Sample"].map(group_map)

    label_mapping = {"mpe": 1, "bpe": 0}
    df_full = df_features[df_features["group"].isin(label_mapping)].copy()
    df_full.sort_values(by="Sample", inplace=True)
    df_full.reset_index(drop=True, inplace=True)

    x_all = df_full[required_features]
    y_all = df_full["group"].map(label_mapping)
    y_pred_prob = pipeline.predict_proba(x_all)[:, 1]
    y_pred_logit = pipeline.decision_function(x_all)
    y_pred_label = (y_pred_prob >= threshold).astype(int)

    tn_value, fp_value, fn_value, tp_value = confusion_matrix(
        y_all, y_pred_label, labels=[0, 1]
    ).ravel()
    auc_value = roc_auc_score(y_all, y_pred_prob)
    print(f"ULP Multi AUC: {auc_value:.4f}")
    print(
        "Confusion Matrix: "
        f"TN={tn_value}, FP={fp_value}, FN={fn_value}, TP={tp_value}"
    )

    result_df = df_full[["Sample", "group"]].copy()
    result_df["Multi_Score_Prob"] = y_pred_prob
    result_df["Multi_Score_Logit"] = y_pred_logit
    result_df["Multi_Predicted_Label"] = y_pred_label
    output_columns = [
        "Sample",
        "Multi_Score_Prob",
        "Multi_Score_Logit",
        "Multi_Predicted_Label",
    ]
    result_df[output_columns].to_csv(PREDICTION_CSV, index=False)


if __name__ == "__main__":
    main()
