# Recorded software versions — author verification pending

The following section is copied from the source working README. It records supplied historical versions, not versions independently verified by executing all pipelines. Resolve discrepancies using analysis logs before publication. WGBS versions still need the actual server environment.

## Software dependencies and versions

The following software versions were used for the manuscript analyses.

### Python

| Software/package | Version |
|---|---:|
| Python | 3.12.9 |
| NumPy | 2.2.4 |
| pandas | 2.2.3 |
| SciPy | 1.15.2 |
| statsmodels | 0.14.4 |
| scikit-learn | 1.6.1 |
| joblib | 1.4.2 |
| matplotlib | 3.10.1 |
| seaborn | 0.13.2 |
| pysam | 0.23.0 |
| tqdm | 4.67.1 |

### Command-line bioinformatics software

| Software | Version |
|---|---:|
| Dorado | 1.3.1 |
| NanoFilt | 2.8.0 |
| seqtk | 1.5-r133 |
| minimap2 | 2.1 |
| samtools | 1.21 |
| mosdepth | 0.3.10 |
| MultiQC | 1.27.1 |
| modkit | 0.4.2 |
| NanoPlot | 1.43.0 |
| pigz | 2.6 |
| bgzip/HTSlib | 1.21 |
| gzip | 1.12 |

### ichorCNA/R environment

| Software/component | Version |
|---|---:|
| R | 4.4.2 |
| ichorCNA | 0.3.2 |
| HMMcopy | 1.48.0 |
| Platform | x86_64-conda-linux-gnu |
| Operating system | CentOS Linux 7 |

### NGS and EPIC v2 preprocessing dependencies

The supplied upstream processing notes report the following historical versions;
they were not verified by executing the complete pipelines in this integration.

| Component | Version / provenance status |
|---|---|
| fastp | 0.19.5 |
| Sentieon | 202010 |
| BWA (provider workflow) | 0.7.17-r1188 |
| samtools (provider workflow) | 1.9; distinct from the Nanopore environment |
| MethylDackel | 0.5.0, as reported in the supplied preprocessing notes |
| gctk | Source and version require author confirmation |
| R (NGS/CHIP preprocessing notes) | 4.5.0; separate from the ichorCNA environment |
| methylKit (5 kb summarization) | 1.34.0; provider notes separately mention 1.16.1 |
| dplyr | 1.1.4 |
| minfi | 1.54.1 |
| IlluminaHumanMethylationEPICv2manifest | 0.3.0 |
| ChAMP / ChAMPdata | Supplied notes specify YuanTian1991 GitHub forks; commits not supplied |
| conumee | Custom EPICv2-capable archive mentioned in source notes; archive/version not supplied |
| GenomicRanges / IRanges / sva | Match the R/Bioconductor environment; exact interval-package versions not supplied; notes report sva 3.56.0 |

Install R/Bioconductor dependencies separately from `requirements.txt`.
The supplied installation sources for ChAMP are `YuanTian1991/ChAMPdata`
and `YuanTian1991/ChAMP`; record exact commits before claiming a pinned environment.
Analysis scripts do not install or update packages automatically.

