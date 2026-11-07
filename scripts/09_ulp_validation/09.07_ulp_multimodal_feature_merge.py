"""Purpose: Merge ULP METH, CNA, and FRAG feature tables.

Edit the paths below, then run this script directly.
"""

import pandas as pd


FRAG_FILE = "results/09_ulp_validation/ULP_FRAG_features.csv"
CNA_FILE = "results/09_ulp_validation/ULP_CNA_features.csv"
METH_FILE = "results/09_ulp_validation/ULP_METH_features.csv"
GROUP_FILE = "metadata/ULP_groupfile.csv"
OUTPUT_FILE = "results/09_ulp_validation/ULP_Multi_features.csv"


def main() -> None:
    df_frag = pd.read_csv(FRAG_FILE)
    df_cna = pd.read_csv(CNA_FILE)
    df_meth = pd.read_csv(METH_FILE)
    group_df = pd.read_csv(GROUP_FILE)

    for df in [df_frag, df_cna, df_meth, group_df]:
        df["Sample"] = df["Sample"].astype(str).str.strip()

    # FRAG-first column order matches the supplied ULP_Multi_features.csv table.
    df_merge = pd.merge(df_frag, df_cna, on="Sample", how="inner")
    df_merge = pd.merge(df_merge, df_meth, on="Sample", how="inner")

    metadata_columns = ["Sample", "group", "split"]
    df_full = pd.merge(df_merge, group_df[metadata_columns], on="Sample", how="inner")
    df_full = df_full[df_full["group"].isin(["mpe", "bpe"])].copy()
    df_full.sort_values(by="Sample", inplace=True)
    df_full.reset_index(drop=True, inplace=True)

    excluded_columns = ["Sample", "group", "split", "subgroup"]
    feature_columns = [column for column in df_full.columns if column not in excluded_columns]
    df_final = df_full[["Sample"] + feature_columns]
    df_final.to_csv(OUTPUT_FILE, index=False)


if __name__ == "__main__":
    main()
