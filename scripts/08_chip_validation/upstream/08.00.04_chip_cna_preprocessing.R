# Purpose: prepare control-normalized EPIC v2 MethylSets for reference and validation arrays.
# Inputs: sample sheets (Basename) and ArrayId-indexed metadata (Sample, Group, Batch).
# Outputs: separate reference and validation MethylSet RDS files.
# Usage: Rscript scripts/08_chip_validation/upstream/08.00.04_chip_cna_preprocessing.R
library(minfi)

REFERENCE_SHEET <- "metadata/chip_preprocessing/reference_sample_sheet.csv"
REFERENCE_METADATA <- "metadata/chip_preprocessing/reference_array_metadata.tsv"
VALIDATION_SHEET <- "metadata/chip_preprocessing/validation_sample_sheet.csv"
VALIDATION_METADATA <- "metadata/chip_preprocessing/validation_array_metadata.tsv"
REFERENCE_RDS <- "results/08_chip_validation/preprocessing/cna/reference_methylset.rds"
VALIDATION_RDS <- "results/08_chip_validation/preprocessing/cna/validation_methylset.rds"

prepare_methylset <- function(sample_sheet_file, meta_file, output_file) {

  meta = read.table(meta_file, header=T, sep="\t",row.names=1, stringsAsFactors=FALSE)
  ## iDat
  sheet = read.csv(sample_sheet_file, stringsAsFactors = FALSE)

  ## RGChannelSet
  RGset <- read.metharray.exp(targets = sheet)

  pData(RGset)$Label = meta[sampleNames(RGset), "Sample"]
  pData(RGset)$Group = meta[sampleNames(RGset), "Group"]
  pData(RGset)$Batch = meta[sampleNames(RGset), "Batch"]
  colnames(RGset) = pData(RGset)$Label


  ## MethylSet & normalization
  Mset = preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls")

  saveRDS(Mset, file=output_file)
}

prepare_methylset(REFERENCE_SHEET, REFERENCE_METADATA, REFERENCE_RDS)
prepare_methylset(VALIDATION_SHEET, VALIDATION_METADATA, VALIDATION_RDS)
