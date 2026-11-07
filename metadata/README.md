# Metadata included in this candidate

Root PE CSVs and CHIP *.example.* tables are synthetic schema examples. They do
not match the real cohort's sample identifiers, reproduce clinical distributions,
or recover manuscript results. Approved real metadata must match the input
features/predictions when running cohort analyses. Sample identifiers must be
unique. group=mpe/bpe is encoded 1/0; preserve the manuscript covariate coding.

Sample/group/split support model evaluation; age/gender support the discovery
GLM; LBC/CB/Combined_Cyto support cytology analyses. The example cytology values
illustrate schema only; confirm actual script-specific positive coding before use.

public_tissue_wgbs/ contains supplied public accessions and provenance manifests,
not synthetic PE records. Tumor/non-malignant tissue groups are not mpe/bpe.
No original participant identifiers or linkage tables are included.
