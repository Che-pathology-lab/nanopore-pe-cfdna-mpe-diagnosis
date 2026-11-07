"""Purpose: Vif filter for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant

INPUT_FILE = "results/05_nanopore_multimodal/discovery_Multi_features_wide_table.csv"
GROUP_FILE = "metadata/discovery_groupfile.csv"
OUTPUT_CLEAN_CSV = "results/05_nanopore_multimodal/discovery_Multi_features_vif_filtered.csv"
VIF_THRESHOLD = 10.0

def calculate_vif_and_drop(X, threshold=10.0):
    X_calc = add_constant(X)
    dropped_features = []
    
    while True:
        vif = pd.Series(
            [variance_inflation_factor(X_calc.values, i) for i in range(X_calc.shape[1])],
            index=X_calc.columns
        ).drop('const')
        
        max_vif = vif.max()
        if max_vif > threshold:
            max_feat = vif.idxmax()
            print(f"delete: {max_feat} (VIF: {max_vif:.2f})")
            X_calc = X_calc.drop(columns=[max_feat])
            dropped_features.append(max_feat)
        else:
            break
            
    return dropped_features

def main():

    df = pd.read_csv(INPUT_FILE)
    group_df = pd.read_csv(GROUP_FILE)
    df = pd.merge(df, group_df[['Sample', 'split']], on='Sample', how='inner').set_index('Sample')
    numeric_df_all = df.select_dtypes(include=[np.number]).replace([np.inf, -np.inf], np.nan).dropna(axis=1)
    X_train = numeric_df_all[df['split'] == 'train']
    to_drop = calculate_vif_and_drop(X_train, VIF_THRESHOLD)
    df_reduced_all = numeric_df_all.drop(columns=to_drop)
    df_reduced_all.to_csv(OUTPUT_CLEAN_CSV)


if __name__ == "__main__":
    main()
