"""Purpose: Extract NGS fragment end-motif features.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np
from scipy.stats import entropy

MOTIF_FILE = "results/07_ngs_validation/FRAG/NGS_endmotif_matrix.csv"
OUTPUT_CSV = "results/07_ngs_validation/FRAG/NGS_FRAG_features_wide_table.csv"

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

def main():
    df_features = extract_endmotif_features(MOTIF_FILE)
    df_features.reset_index(inplace=True)
    df_features.to_csv(OUTPUT_CSV, index=False)


if __name__ == "__main__":
    main()
