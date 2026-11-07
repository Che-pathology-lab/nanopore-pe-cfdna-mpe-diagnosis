"""Purpose: Calculate ULP FRAG features using the discovery feature definitions.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np
from scipy.stats import entropy

MOTIF_FILE = "results/09_ulp_validation/ULP_endmotif_matrix.csv"
LENGTH_FILE = "results/09_ulp_validation/ULP_lengths_matrix.csv"
OUTPUT_CSV = "results/09_ulp_validation/ULP_FRAG_features.csv"
LENGTH_RANGES = [(220, 270)]

def extract_endmotif_features(filepath):

    df_raw = pd.read_csv(filepath, index_col=0)
    df_raw.index.name = "Sample"
    df_raw.index = df_raw.index.astype(str)
    motif_4mer_cols = [c for c in df_raw.columns if len(c) == 4 and set(c).issubset(set('ACGT'))]
    features_dict = {}
    for col in motif_4mer_cols:
        features_dict[col] = df_raw[col]
    bases = ['A', 'C', 'G', 'T']

    trimers = [b1+b2+b3 for b1 in bases for b2 in bases for b3 in bases]
    for trimer in trimers:
        cols = [c for c in motif_4mer_cols if c.startswith(trimer)]
        if cols: features_dict[trimer] = df_raw[cols].sum(axis=1)

    dimers = [b1+b2 for b1 in bases for b2 in bases]
    for dimer in dimers:
        cols = [c for c in motif_4mer_cols if c.startswith(dimer)]
        if cols: features_dict[dimer] = df_raw[cols].sum(axis=1)

    for monomer in bases:
        cols = [c for c in motif_4mer_cols if c.startswith(monomer)]
        if cols: features_dict[monomer] = df_raw[cols].sum(axis=1)

    entropy_list = []
    for _, row in df_raw.iterrows():
        counts = row[motif_4mer_cols].values.astype(float)
        if counts.sum() > 0:
            probs = counts / counts.sum()
            ent = entropy(probs, base=2)
        else:
            ent = np.nan
        entropy_list.append(ent)
        
    features_dict['Motif_Entropy'] = entropy_list

    df_features = pd.DataFrame(features_dict, index=df_raw.index)
    return df_features

def extract_length_features(filepath, ranges):

    df_raw = pd.read_csv(filepath, index_col=0)
    df_raw.index.name = "Sample"
    df_raw.index = df_raw.index.astype(str)
    all_cols = df_raw.columns
    
    features_dict = {}
    
    for target_start, target_end in ranges:
        target_cols = [c for c in all_cols if c.isdigit() and target_start <= int(c) <= target_end]
        feature_name = f'Freq_{target_start}_{target_end}'
        features_dict[feature_name] = df_raw[target_cols].sum(axis=1)
    df_features = pd.DataFrame(features_dict, index=df_raw.index)
    return df_features

def main():
    df_motif = extract_endmotif_features(MOTIF_FILE)
    df_length = extract_length_features(LENGTH_FILE, LENGTH_RANGES)
    df_combined = df_motif.merge(df_length, left_index=True, right_index=True, how='inner')

    out_path = OUTPUT_CSV
    df_combined.reset_index(inplace=True)
    df_combined.to_csv(out_path, index=False)


if __name__ == "__main__":
    main()
