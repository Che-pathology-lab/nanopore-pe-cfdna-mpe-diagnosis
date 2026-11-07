# Purpose: average normalized EPIC v2 probe beta values within supplied Nanopore DMRs.
# Inputs: normalized ChAMP workspace, hg38 probe annotation, 0-based DMR intervals.
# Output: 0-1 DMR means; the original overlap and mean rules are retained.
# Usage: Rscript scripts/08_chip_validation/upstream/08.00.02_chip_dmr_summary.R
library(dplyr)

ANNOTATION_FILE <- "references/epic_v2_hg38_probe_positions.tsv"
NORMALIZED_FILE <- "results/08_chip_validation/preprocessing/methylation/step2-champ_myNorm.Rdata"
DMR_FILE <- "references/chip_nanopore_dmr_intervals.tsv"
OUTPUT_FILE <- "results/08_chip_validation/preprocessing/methylation/chip_dmr_beta.tsv"

gplhg38 <- read.delim(ANNOTATION_FILE)

load(file = NORMALIZED_FILE)
mycombat <- as.data.frame(myCombat)
mycombat$id_cpg <- rownames(mycombat)
mycombat <- mycombat[, c("id_cpg", setdiff(names(mycombat), "id_cpg"))]
rownames(mycombat) <- NULL

result <- inner_join(
  x = gplhg38,
  y = mycombat,
  by = "id_cpg"
) %>% 

  arrange(match(id_cpg, gplhg38$id_cpg))

train_combined <- read.delim(DMR_FILE)
train_combined$ID <- paste0(train_combined$chr, "_",train_combined$start)
train_combined <- train_combined[, c("ID", setdiff(names(train_combined), "ID"))]
train_combined$start <- train_combined$start + 1



result1_chr <- result %>% filter(chr == "chr1")
result1_chr$chr <- as.character(result1_chr$chr)

result_filtered_chr <- result1_chr %>%

  inner_join(
    train_combined %>% dplyr::select(chr, dmr_start = start, dmr_end = end),
    by = "chr"
  ) %>%

  filter(position >= dmr_start & position <= dmr_end) %>%

  dplyr::select(-dmr_start, -dmr_end) 

result_filtered <- result_filtered_chr
chrs <- paste0("chr", 2:22)
for (i in chrs) {

  result1_chr <- result %>% filter(chr == i)
  result1_chr$chr <- as.character(result1_chr$chr)

  result_filtered_chr <- result1_chr %>%

    inner_join(
      train_combined %>% dplyr::select(chr, dmr_start = start, dmr_end = end),
      by = "chr"
    ) %>%

    filter(position >= dmr_start & position <= dmr_end) %>%

    dplyr::select(-dmr_start, -dmr_end) 
  result_filtered <- rbind(result_filtered, result_filtered_chr)
}


cpg_dmr_joined <- result_filtered %>%

  inner_join(
    train_combined %>% dplyr::select(DMR_ID = ID, dmr_chr = chr, dmr_start = start, dmr_end = end),
    by = c("chr" = "dmr_chr")
  ) %>%

  filter(position >= dmr_start & position <= dmr_end)

cpg_dmr_joined <- cpg_dmr_joined[,-(2:3)]

dmr_mean <- cpg_dmr_joined %>%

  group_by(DMR_ID,dmr_start, dmr_end) %>%

  summarise(
    across(
      .cols = -c(id_cpg),
      .fns = mean,
      .names = "{.col}"
    ),
    .groups = "drop"
  )

colnames(dmr_mean)[[1]] <- "ID"
colnames(dmr_mean)[[2]] <- "start"
colnames(dmr_mean)[[3]] <- "end"

dmr_mean$chr <- sub("_.*", "", dmr_mean$ID)
dmr_mean <- dmr_mean[, c("ID", "chr", setdiff(names(dmr_mean), c("ID", "chr")))]
dmr_mean$start <- dmr_mean$start-1
write.table(dmr_mean, OUTPUT_FILE,
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "NA")
