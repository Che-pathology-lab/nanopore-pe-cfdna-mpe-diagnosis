# Reproducibility scope and open items

This copy provides analysis code and an executable illustrative example. It does
not claim raw-data-to-result reproduction or manuscript numerical revalidation.
Synthetic examples supplement, rather than replace, access to the data required
to evaluate the study. Reviewers may need the real approved inputs and models.

## Software

requirements.txt records Python dependencies; manuscript Python is 3.12.9.
External tools used across modules include Dorado, minimap2, samtools, modkit,
NanoFilt/NanoPlot, mosdepth, MultiQC, seqtk, pigz, bgzip, ichorCNA/HMMcopy;
NGS adds fastp, Sentieon, MethylDackel/gctk and methylKit; EPIC adds ChAMP,
minfi and a custom conumee build. WGBS adds SRA Toolkit, Cutadapt, BSMAP,
Picard/Java, MOABS/mcall and bedtools. Obtain the actual historical versions,
custom sources and reference checksums. A Python-only install covers the example,
not all sequencing pipelines. The source README records some versions but several
upstream builds remain unresolved; version claims should be checked against run logs.

## Scientific and input gaps retained from source

- DMR missingness filtering currently precedes training-sample selection.
  Compare its effect before aligning a train-only Methods description.
- Hyperparameters are selected on the whole training set before OOF predictions;
  feature screening is not fully nested. Do not call this unbiased nested CV.
- Final Multi coefficients have 19 candidate entries and 10 nonzero coefficients.
  This is not the feature count of every separate single-modality model.
- Cluster 2 labels are manually assigned; the script evaluates supplied labels
  rather than reproducing their assignment or external cluster projection.
- NGS historical inputs have a 9-region discrepancy in the supplied overlap check.
- CHIP identity mapping and the historical 40-bin exclusion rule are unresolved;
  a custom conumee build is missing. Include the frozen interval template only
  after provenance/release review and retain an explicit explanation.
- Several ULP upstream matrices, repeated-PoN inputs and provider preprocessing
  commands are not supplied. References to them are not self-contained inputs.
- WGBS is full-depth regional methylation/CNA processing only. Some library
  orientation/QC assumptions and exact environment builds remain unresolved.

Scientific scripts have not been corrected to resolve these issues because doing
so can change results. Resolve them with the authors and describe accurately.
Human data access restrictions require a specific access mechanism, not just
an unspecified 'available on request' sentence. No manuscript-compliance claim
is made until the target journal, data access and central artifact access are settled.
