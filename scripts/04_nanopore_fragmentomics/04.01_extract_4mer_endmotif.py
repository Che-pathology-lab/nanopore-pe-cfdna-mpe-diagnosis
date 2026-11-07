"""Purpose: Extract 4mer endmotif for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pysam
import pandas as pd
import os
import glob
from multiprocessing import Pool


INPUT_DIR = "results/01_nanopore_preprocessing/1.4_mapping"
OUTPUT_FILE = "results/04_nanopore_fragmentomics/4mer_endmotif_matrix.csv"
REF_FASTA_PATH = "references/hg38.fa" 
MIN_MAPQ = 20
MOTIF_LEN = 4
NUM_PROCESSES = 10
COMPLEMENT_MAP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def reverse_complement(seq):
    return seq.translate(COMPLEMENT_MAP)[::-1]

def get_end_motifs(args):
    bam_path, ref_fasta_path = args
    sample_name = os.path.basename(bam_path).replace(".bam", "")
    motif_counts = {}
    
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam, \
         pysam.FastaFile(ref_fasta_path) as ref_genome:
        
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < MIN_MAPQ:
                continue
            
            chrom = read.reference_name
            
            if not read.is_reverse:
                start_pos = read.reference_start
                motif = ref_genome.fetch(chrom, start_pos, start_pos + MOTIF_LEN)
                motif = motif.upper()
                
            else:
                end_pos = read.reference_end
                seq_raw = ref_genome.fetch(chrom, end_pos - MOTIF_LEN, end_pos)
                motif = reverse_complement(seq_raw).upper()
                
            if len(motif) == MOTIF_LEN and "N" not in motif:
                motif_counts[motif] = motif_counts.get(motif, 0) + 1

    return sample_name, motif_counts

def main():
    bam_files = glob.glob(os.path.join(INPUT_DIR, "*.bam"))
    task_args = [(f, REF_FASTA_PATH) for f in bam_files]
    all_data = {}
    
    with Pool(processes=NUM_PROCESSES) as pool:
        results = pool.imap_unordered(get_end_motifs, task_args)
        
        count = 0
        for sample_name, motifs in results:
            count += 1
            print(f"[{count}/{len(bam_files)}] Processed: {sample_name}")
            
            if motifs is not None:
                all_data[sample_name] = motifs

    df = pd.DataFrame.from_dict(all_data, orient='index')
    df = df.fillna(0)
    df = df.reindex(sorted(df.columns), axis=1)
    df_freq = df.div(df.sum(axis=1), axis=0)
    df_freq.to_csv(OUTPUT_FILE, index_label="Sample")


if __name__ == "__main__":
    main()
