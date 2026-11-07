# Purpose: fit EPIC v2 CNA reference profiles and export conumee bin ratios.
# Inputs: control-normalized reference/validation MethylSet RDS files.
# Outputs: per-sample .bins.txt files containing ratio and shift_ratio.
# Usage: Rscript scripts/08_chip_validation/upstream/08.00.05_chip_cna_bins.R
library(conumee)

REFERENCE_RDS <- "results/08_chip_validation/preprocessing/cna/reference_methylset.rds"
VALIDATION_RDS <- "results/08_chip_validation/preprocessing/cna/validation_methylset.rds"
OUTPUT_DIR <- "results/08_chip_validation/preprocessing/cna/bins"

# Retain the supplied EPICv2 annotation call. Confirm genome build, bin size
# and the exact custom conumee build before attempting full reproduction.
ann.cnv = CNV.create_anno(array_type = "EPICv2")

mset.tumor = readRDS(VALIDATION_RDS)
mset.control = readRDS(REFERENCE_RDS)

mset.probes = intersect(rownames(mset.tumor), rownames(mset.control))
common.probes = intersect(mset.probes, names(ann.cnv@probes))

mset.tumor = mset.tumor[common.probes,]
mset.control = mset.control[common.probes,]
ann.cnv@probes = ann.cnv@probes[common.probes]

cat("Common probes: ", length(common.probes), "\n\n")

cnv.tumor = CNV.load(mset.tumor)
cnv.control = CNV.load(mset.control)


for( x in names(cnv.tumor)){ 

  fit = CNV.fit(cnv.tumor[x], cnv.control, ann.cnv)
  fit = CNV.bin(fit)

  bin.res = data.frame(fit@anno@bins, ratio=fit@bin$ratio, shift_ratio=fit@bin$ratio-fit@bin$shift)
  write.table(bin.res, file=file.path(OUTPUT_DIR, paste0(x,".bins.txt")), quote=F, sep="\t", row.names=F)
}
