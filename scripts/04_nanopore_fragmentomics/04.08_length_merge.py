"""Purpose: Length merge for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np
import os
import glob


INPUT_ROOT_DIR = "results/04_nanopore_fragmentomics/nanoplot"
OUTPUT_FILE = "results/04_nanopore_fragmentomics/nanopore_lengths_matrix.csv"
MIN_LEN = 0
MAX_LEN = 3000
DATA_FILE_PATTERN = "*NanoPlot-data.tsv.gz"


def extract_lengths_from_tsv(file_path):
    df = pd.read_csv(file_path, sep='\t', compression='gzip')
    lengths = df['aligned_lengths'].dropna().astype(int).tolist()
    valid_lengths = [l for l in lengths if MIN_LEN <= l <= MAX_LEN]  
    return valid_lengths

def main():
    search_pattern = os.path.join(INPUT_ROOT_DIR, "**", DATA_FILE_PATTERN)
    data_files = glob.glob(search_pattern, recursive=True)
    all_sample_data = {}

    for file_path in data_files:
        filename = os.path.basename(file_path)
        sample_name = filename.replace("_NanoPlot-data.tsv.gz", "")
        
        lengths = extract_lengths_from_tsv(file_path)
        counts, _ = np.histogram(lengths, bins=range(MIN_LEN, MAX_LEN + 2))
        all_sample_data[sample_name] = counts
    df = pd.DataFrame.from_dict(all_sample_data, orient='index')
    df.columns = [str(i) for i in range(MIN_LEN, MAX_LEN + 1)]
    df = df.sort_index()
    df_freq = df.div(df.sum(axis=1), axis=0)
    df_freq.to_csv(OUTPUT_FILE, index_label="Sample")

if __name__ == "__main__":
    main()
