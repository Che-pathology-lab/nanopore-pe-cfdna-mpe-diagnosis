# Purpose: tile GM-seq CpG counts into 5 kb bins and calculate methylation/CpG counts.
# Inputs: seven-column CpG tables and metadata/ngs_groupfile.csv (Sample, group).
# Outputs: count and 0-1 methylation tables in preprocessing/methylation.
# Usage: Rscript scripts/07_ngs_validation/upstream/07.00.05_methylkit_5kb.R
library(methylKit)
library(GenomicRanges)
library(IRanges)

INPUT_DIR <- "results/07_ngs_validation/preprocessing/cpg_methylkit"
GROUP_FILE <- "metadata/ngs_groupfile.csv"
OUTPUT_DIR <- "results/07_ngs_validation/preprocessing/methylation"


sample_info <- read.csv(GROUP_FILE, check.names = FALSE)
sample_info$ID <- sample_info$Sample
rownames(sample_info) <- NULL




file_paths <- file.path(INPUT_DIR, paste0(sample_info$ID, ".CpG.txt"))
file.list <- as.list(file_paths)

group=ifelse(sample_info$group == "mpe", 1, 0)
id.list=as.list(sample_info$ID)

myobj <- methRead(
  file.list,
  sample.id = id.list,
  treatment = group,
  context = "CpG",
  assembly = "hg38",
  sep = "\t",
  header = TRUE,
  mincov = 1
)
regions <- tileMethylCounts(myobj, win.size = 5000, step.size = 5000)
meth2 <- methylKit::unite(regions, destrand=FALSE, min.per.group=1L)

total_needed <- length(id.list) * 3
new_col_names <- unlist(lapply(id.list, function(id) {
  c(paste0(id, "_coverage"), 
    paste0(id, "_numCs"), 
    paste0(id, "_numTs"))
}))
colnames(meth2)[5:(4 + total_needed)] <- new_col_names
meth2_df <- as.data.frame(meth2)

# All outputs use repository-relative paths; no working-directory changes.
write.table(meth2_df, file.path(OUTPUT_DIR, "meth_region5k.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            na = "NA")
################################################
meth2_df <- read.delim(file.path(OUTPUT_DIR, "meth_region5k.txt"), 
                       sep = "\t", 
                       header = TRUE, 
                       check.names = FALSE,
                       stringsAsFactors = FALSE)


numCs_cols <- grep(pattern = "_numCs$", x = colnames(meth2_df), value = TRUE)

prefixes <- sub(pattern = "_numCs$", replacement = "", x = numCs_cols)
coverage_cols <- paste0(prefixes, "_coverage")

valid_idx <- coverage_cols %in% colnames(meth2_df)
prefixes_valid <- prefixes[valid_idx]
numCs_cols_valid <- numCs_cols[valid_idx]
coverage_cols_valid <- coverage_cols[valid_idx]

meth_list <- lapply(seq_along(prefixes_valid), function(i) {
  numC <- meth2_df[[numCs_cols_valid[i]]]
  cov <- meth2_df[[coverage_cols_valid[i]]]

  meth_rate <- round(ifelse(cov == 0, NA, numC / cov), digits = 2)
  return(meth_rate)
})


names(meth_list) <- paste0(prefixes_valid, "_meth")
meth2_df <- cbind(meth2_df, do.call(cbind, meth_list))

head(meth2_df)
pattern <- "_numCs|_numTs"
cols_to_delete <- grep(pattern = pattern, x = colnames(meth2_df), value = TRUE)

if (length(cols_to_delete) == 0) {
  warning("No _numCs or _numTs columns found; output columns are retained.")
} else {

  meth2_df <- meth2_df[, !colnames(meth2_df) %in% cols_to_delete, drop = FALSE]
  cat("Removed count columns after calculating methylation rates.\n")
}
meth2_df$ID <- paste(meth2_df$chr, meth2_df$start, sep = "_")

meth2_df2 <- meth2_df
meth2_df2$ID <- paste0("chr", meth2_df2$ID)
meth2_df2 <- meth2_df2[, c("ID", setdiff(names(meth2_df2), "ID"))]

write.table(meth2_df2, file.path(OUTPUT_DIR, "meth_region5k2.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            na = "NA")




sample_names <- prefixes



cpg_gr_list <- lapply(myobj, function(sample_data) {

  cpg_df <- as.data.frame(sample_data)

  GRanges(
    seqnames = cpg_df$chr,
    ranges = IRanges(start = cpg_df$start, end = cpg_df$end),
    strand = "*"
  )
})


names(cpg_gr_list) <- sample_names

meth_gr <- GRanges(
  seqnames = meth2_df2$chr,
  ranges = IRanges(start = meth2_df2$start, end = meth2_df2$end),
  strand = "*"
)


for (sample in sample_names) {

  sample_cpg_gr <- cpg_gr_list[[sample]]


  dmr_cpg_count <- countOverlaps(
    query = meth_gr,
    subject = sample_cpg_gr
  )


  meth2_df2[[paste0(sample, "_cpg")]] <- dmr_cpg_count
}

write.table(meth2_df2, file.path(OUTPUT_DIR, "hg38_region5k.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "NA")
