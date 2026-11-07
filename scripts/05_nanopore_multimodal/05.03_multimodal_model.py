"""Purpose: Multimodal model for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
import joblib

INPUT_FILE = "results/05_nanopore_multimodal/discovery_Multi_features_vif_filtered.csv"
GROUP_FILE = "metadata/discovery_groupfile.csv"
MODEL_OUTPUT_PATH = "results/05_nanopore_multimodal/discovery_Multi_logistic_model.joblib"
PREDICTION_CSV = "results/05_nanopore_multimodal/discovery_Multi_prediction_results_with_scores.csv"
COEF_CSV = "results/05_nanopore_multimodal/discovery_Multi_coefficients.csv"
PRIMARY_RANDOM_STATE = 100

def main():

    df_features = pd.read_csv(INPUT_FILE)
    df_split = pd.read_csv(GROUP_FILE)
    df_full = pd.merge(df_features, df_split[['Sample', 'group', 'split']], on='Sample', how='inner')
    feature_cols = [col for col in df_features.columns if col != 'Sample']
    df_full = df_full[df_full['group'].isin(['mpe', 'bpe'])].copy()
    df_full.sort_values(by='Sample', inplace=True)
    df_full.reset_index(drop=True, inplace=True)
    X_all = df_full[feature_cols]
    y_all = df_full['group'].map({'mpe': 1, 'bpe': 0})
    train_mask = df_full['split'] == 'train'
    test_mask = df_full['split'] == 'test'

    X_train = X_all[train_mask].copy()
    y_train = y_all[train_mask].copy()
    X_test = X_all[test_mask].copy()
    y_test = y_all[test_mask].copy()
    
    train_df = df_full[train_mask].copy()
    test_df = df_full[test_mask].copy()

    pipe = Pipeline([
        ('scaler', StandardScaler()), 
        ('clf', LogisticRegression(
            solver='liblinear', random_state=PRIMARY_RANDOM_STATE,
            class_weight='balanced', max_iter=5000
        ))
    ])

    param_grid = {
        'clf__C': [0.001, 0.01, 0.1, 1, 10, 100],
        'clf__penalty': ['l1','l2']
    }
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=PRIMARY_RANDOM_STATE)
    
    print("\nStarting GridSearchCV...")
    grid_search = GridSearchCV(pipe, param_grid, cv=cv, scoring='roc_auc', n_jobs=-1, refit=True)
    grid_search.fit(X_train, y_train)
    best_model = grid_search.best_estimator_
    print(f"Best CV AUC: {grid_search.best_score_:.4f}")
    print(f"\nBest params: {grid_search.best_params_}")

    y_train_oof_prob = cross_val_predict(
        best_model, X_train, y_train, cv=cv, method='predict_proba', n_jobs=-1
    )[:, 1]

    y_train_oof_logit = cross_val_predict(
        best_model, X_train, y_train, cv=cv, method='decision_function', n_jobs=-1
    )

    auc_train_oof = roc_auc_score(y_train, y_train_oof_prob)
    fpr_train, tpr_train, thresholds_train = roc_curve(y_train, y_train_oof_prob)
    optimal_idx = np.argmax(tpr_train - fpr_train)
    optimal_threshold = thresholds_train[optimal_idx]
    
    print(f"Train OOF AUC: {auc_train_oof:.4f}")
    print(f"Optimal Threshold: {optimal_threshold:.4f}")
    model_package = {
        'pipeline': best_model,
        'threshold': optimal_threshold,
        'features': feature_cols,
        'description': 'Multi Logistic Regression trained'
    }
    joblib.dump(model_package, MODEL_OUTPUT_PATH)

    y_test_pred_prob = best_model.predict_proba(X_test)[:, 1]
    y_test_pred_logit = best_model.decision_function(X_test)
    y_test_pred_label = (y_test_pred_prob >= optimal_threshold).astype(int)

    tn, fp, fn, tp = confusion_matrix(y_test, y_test_pred_label).ravel()
    auc_test = roc_auc_score(y_test, y_test_pred_prob)
    print(f"Test AUC: {auc_test:.4f}")
    print(f"Test Confusion Matrix: TN={tn}, FP={fp}, FN={fn}, TP={tp}")
    train_df = train_df.copy()
    train_df['Multi_Score_Prob'] = y_train_oof_prob
    train_df['Multi_Score_Logit'] = y_train_oof_logit
    train_df['Multi_Predicted_Label'] = (y_train_oof_prob >= optimal_threshold).astype(int)
    train_df['Set'] = 'Train'

    test_df = test_df.copy()
    test_df['Multi_Score_Prob'] = y_test_pred_prob
    test_df['Multi_Score_Logit'] = y_test_pred_logit
    test_df['Multi_Predicted_Label'] = y_test_pred_label
    test_df['Set'] = 'Test'

    final_res = pd.concat([train_df, test_df], ignore_index=True)
    cols = ['Sample', 'Multi_Score_Prob', 'Multi_Score_Logit', 'Multi_Predicted_Label']
    final_res[cols].to_csv(PREDICTION_CSV, index=False)

    clf = best_model.named_steps['clf']
    coefs = clf.coef_[0]
    coef_df = pd.DataFrame({'Feature': feature_cols, 'Coefficient': coefs, 'Abs_Coefficient': abs(coefs)})
    coef_df = coef_df.sort_values(by='Abs_Coefficient', ascending=False)
    coef_df.to_csv(COEF_CSV, index=False)


if __name__ == "__main__":
    main()
