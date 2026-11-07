"""Purpose: Methylation model for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
import joblib


INPUT_FILE = "results/02_nanopore_methylation/discovery_METH_features_vif_filtered.csv"
GROUP_FILE = "metadata/discovery_groupfile.csv"
MODEL_OUTPUT_PATH = "results/02_nanopore_methylation/discovery_METH_logistic_model.joblib"
FEATURE_COLS = ['Global_Hyper_Mean', 'Global_Hypo_Mean']
OUTPUT_RESULT = "results/02_nanopore_methylation/discovery_METH_prediction_results_with_scores.csv"
PRIMARY_RANDOM_STATE = 100

def main():
    df_matrix = pd.read_csv(INPUT_FILE)
    df_split = pd.read_csv(GROUP_FILE)
    df_full = pd.merge(df_matrix, df_split[['Sample', 'group', 'split']], on='Sample', how='inner')
    df_full.sort_values(by='Sample', inplace=True)
    df_full.reset_index(drop=True, inplace=True)
    X_all = df_full[FEATURE_COLS]
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
            solver='liblinear', 
            random_state=PRIMARY_RANDOM_STATE,
            class_weight='balanced', 
            max_iter=5000
        ))
    ])

    param_grid = {
        'clf__C': [0.001, 0.01, 0.1, 1, 10, 100],
        'clf__penalty': ['l1','l2']
    }

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=PRIMARY_RANDOM_STATE)

    grid_search = GridSearchCV(
        estimator=pipe, 
        param_grid=param_grid, 
        cv=cv, 
        scoring='roc_auc', 
        n_jobs=-1, 
        refit=True
    )
    grid_search.fit(X_train, y_train)
    
    best_model = grid_search.best_estimator_

    y_train_oof_prob = cross_val_predict(
        best_model, X_train, y_train, cv=cv, method='predict_proba', n_jobs=-1
    )[:, 1]

    y_train_oof_logit = cross_val_predict(
        best_model, X_train, y_train, cv=cv, method='decision_function', n_jobs=-1
    )

    auc_train_oof = roc_auc_score(y_train, y_train_oof_prob)
    print(f"Train OOF AUC: {auc_train_oof:.4f}")
    fpr_train, tpr_train, thresholds_train = roc_curve(y_train, y_train_oof_prob)
    youden_index = tpr_train - fpr_train
    optimal_idx = youden_index.argmax()
    optimal_threshold = thresholds_train[optimal_idx]
    print(f"Optimal Threshold (from Train OOF): {optimal_threshold:.4f}")

    model_package = {
        'pipeline': best_model,
        'threshold': optimal_threshold,
        'features': FEATURE_COLS,
        'description': 'METH Logistic Regression trained'
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
    train_df['METH_Score_Prob'] = y_train_oof_prob
    train_df['METH_Score_Logit'] = y_train_oof_logit
    train_df['METH_Predicted_Label'] = (y_train_oof_prob >= optimal_threshold).astype(int)
    train_df['Set'] = 'Train'

    test_df = test_df.copy()
    test_df['METH_Score_Prob'] = y_test_pred_prob
    test_df['METH_Score_Logit'] = y_test_pred_logit
    test_df['METH_Predicted_Label'] = y_test_pred_label
    test_df['Set'] = 'Test'

    final_res = pd.concat([train_df, test_df], ignore_index=True)
    

    cols = ['Sample', 'METH_Score_Prob', 'METH_Score_Logit', 'METH_Predicted_Label']

    final_res[cols].to_csv(OUTPUT_RESULT, index=False)

    clf = best_model.named_steps['clf']
    coefs = clf.coef_[0]
    print("\nModel Coefficients:")
    for name, val in zip(FEATURE_COLS, coefs):
        print(f"{name:<20}: {val:.4f}")


if __name__ == "__main__":
    main()
