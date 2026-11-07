"""Purpose: Cna sigtest for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


FEATURE_FILE = "results/03_nanopore_cna/discovery_CNA_features_wide_table.csv"
GROUP_FILE = "metadata/discovery_groupfile.csv"
STATS_FILE = "results/03_nanopore_cna/discovery_CNA_feature_statistics.csv"  
FILTERED_MATRIX_FILE = "results/03_nanopore_cna/discovery_CNA_features_significant.csv"

def main():

    FDR_THRESHOLD = 0.01

    df_features = pd.read_csv(FEATURE_FILE)
    df_split = pd.read_csv(GROUP_FILE)
    df_full = pd.merge(df_features, df_split[['Sample', 'group', 'split']], on='Sample', how='inner')


    exclude_cols = ['Sample', 'group', 'split']
    feature_cols = [c for c in df_full.columns if c not in exclude_cols]
    train_df = df_full[df_full['split'] == 'train']
    mpe_df = train_df[train_df['group'] == 'mpe']
    bpe_df = train_df[train_df['group'] == 'bpe']
    
    results = []
    
    for feat in feature_cols:
        c_vals = mpe_df[feat].values
        b_vals = bpe_df[feat].values
    
        c_median = np.median(c_vals)
        b_median = np.median(b_vals)
    
        combined_vals = np.concatenate([c_vals, b_vals])
        
        if len(np.unique(combined_vals)) <= 1:
            p_val = 1.0
        else:
            try:
                stat, p_val = mannwhitneyu(c_vals, b_vals, alternative='two-sided')
            except ValueError:
                p_val = 1.0
        
        results.append({
            'Feature': feat,
            'P_Value': p_val,
            'MPE_median': c_median,
            'BPE_median': b_median,
            'Diff': c_median - b_median
        })
        
    stats_df = pd.DataFrame(results)

    _, p_adj, _, _ = multipletests(stats_df['P_Value'], method='fdr_bh')
    stats_df['FDR'] = p_adj
    stats_df['Significant'] = stats_df['FDR'] < FDR_THRESHOLD
    stats_df = stats_df.sort_values(by='P_Value')
    stats_df.to_csv(STATS_FILE, index=False)

    sig_features = stats_df[stats_df['Significant']]['Feature'].tolist()
    final_cols = ['Sample'] + sig_features
    df_filtered = df_features[final_cols]
    df_filtered.to_csv(FILTERED_MATRIX_FILE, index=False)

if __name__ == "__main__":
    main()
