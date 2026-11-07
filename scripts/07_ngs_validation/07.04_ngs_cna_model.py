"""Purpose: Ngs cna model for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import joblib
from sklearn.metrics import roc_auc_score, confusion_matrix

MODEL_FILE = "results/03_nanopore_cna/discovery_CNA_logistic_model.joblib"
FEATURE_FILE = "results/07_ngs_validation/CNA/NGS_CNA_features_wide_table.csv"
GROUP_FILE = "metadata/ngs_groupfile.csv" 
PREDICTION_CSV = "results/07_ngs_validation/CNA/NGS_CNA_prediction_results_with_scores.csv"

def main():


    model_pkg = joblib.load(MODEL_FILE)
    pipeline = model_pkg['pipeline']
    threshold = model_pkg['threshold']
    required_features = model_pkg['features']
    df_features = pd.read_csv(FEATURE_FILE)
    df_split = pd.read_csv(GROUP_FILE)

    group_map = dict(zip(df_split['Sample'], df_split['group']))
    df_features['group'] = df_features['Sample'].map(group_map)

    label_mapping = {'mpe': 1, 'bpe': 0}
    df_full = df_features.copy()
    df_full = df_full[df_full['group'].isin(label_mapping.keys())].copy()
    df_full.sort_values(by='Sample', inplace=True)
    df_full.reset_index(drop=True, inplace=True)
    X_all = df_full[required_features]
    y_all = df_full['group'].map(label_mapping)
    y_pred_prob = pipeline.predict_proba(X_all)[:, 1]
    y_pred_logit = pipeline.decision_function(X_all)

    y_pred_label = (y_pred_prob >= threshold).astype(int)
    cm = confusion_matrix(y_all, y_pred_label, labels=[0, 1])
    tn_te, fp_te, fn_te, tp_te = cm.ravel()

    auc_test = roc_auc_score(y_all, y_pred_prob)
    print(f"NGS AUC: {auc_test:.4f}")
    print(f"Confusion Matrix: TN={tn_te}, FP={fp_te}, FN={fn_te}, TP={tp_te}")
    res_df = df_full[['Sample', 'group']].copy()
    
    res_df['CNA_Score_Prob'] = y_pred_prob
    res_df['CNA_Score_Logit'] = y_pred_logit
    res_df['CNA_Predicted_Label'] = y_pred_label

    columns_to_save = ['Sample', 'CNA_Score_Prob', 'CNA_Score_Logit', 'CNA_Predicted_Label']
    res_df[columns_to_save].to_csv(PREDICTION_CSV, index=False)


if __name__ == "__main__":
    main()
