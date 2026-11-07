# Purpose: load EPIC v2 IDATs, perform ChAMP QC, BMIQ and Slide ComBat.
# Inputs: local IDAT directory with Samplesheet.csv (see module README).
# Outputs: ChAMP workspaces; preprocessing output is intended for local use only.
# Usage: Rscript scripts/08_chip_validation/upstream/08.00.01_chip_methylation_preprocessing.R
library(ChAMPdata)
data(AnnoEPICv2)
library(ChAMP)
library(IlluminaHumanMethylationEPICv2manifest)

IDAT_DIR <- "data/chip/idat"
OUTPUT_DIR <- "results/08_chip_validation/preprocessing/methylation"
CORES <- 8

myLoad <- champ.load(directory = IDAT_DIR, arraytype = "EPICv2")

save(myLoad, file = file.path(OUTPUT_DIR, "step1-output.Rdata"))
champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Group)

myNorm <- champ.norm(beta = myLoad$beta, arraytype="EPICv2", cores = CORES, method = "BMIQ")

myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname="Slide")

dim(myCombat) 
pD=myLoad$pd

save(myNorm,myCombat,pD,file = file.path(OUTPUT_DIR, "step2-champ_myNorm.Rdata"))
