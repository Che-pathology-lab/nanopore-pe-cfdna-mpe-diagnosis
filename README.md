# CURE multimodal PE cfDNA analysis

Analysis code associated with the Nanopore multimodal framework for malignant
pleural effusion diagnosis. The repository contains preprocessing, feature
extraction, model development/evaluation and figure/metric calculations.

## Release scope

This release includes analysis code, recorded dependencies, synthetic PE metadata
schemas, supplied public tissue WGBS accession manifests and a runnable synthetic
ROC/bootstrap example. It does not include patient-level clinical records,
participant features/predictions, cluster assignments, trained model packages,
raw sequencing/IDAT files or reference data with unresolved release provenance.
The original working repository and its study data have not been modified.

Synthetic examples demonstrate execution, not the manuscript results. Reproducing
the cohort results requires approved study data, exact references and model
artifacts. Required inputs and known limitations are documented in
[reproducibility scope](docs/reproducibility_scope.md),
[software versions](docs/software_versions.md) and [references](references/README.md).

## Repository structure

- scripts/: numbered analysis scripts, with editable paths/settings.
- metadata/: synthetic PE/array schemas and public tissue WGBS accessions.
- examples/: invented matching metadata/predictions and a runnable example.
- docs/: software records, execution scope and verification notes.
- references/: instructions for obtaining the exact required references.

## Installation and synthetic example

The manuscript Python version is 3.12.9. In an isolated compatible environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python examples/run_synthetic_roc_ci.py
```

The example uses 20 invented samples and reuses modules 11 and 15. It writes
examples/output/synthetic_roc.pdf and synthetic_metrics.csv (24 metric rows).
Expected rounded metric output is examples/expected_synthetic_metrics.csv.
The demo uses 200 bootstrap iterations for speed; original analysis scripts
retain their settings (5000 in module 15, random seed 2026). No fitted model
or real cohort data are loaded. The generated figure is explicitly synthetic.

## Analysis workflow

Run manuscript scripts from the repository root after editing paths at the top.
The original scientific commands, model algorithms and thresholds are retained.

| Modules | Function |
|---|---|
| 01–05 | ONT preprocessing, METH/CNA/FRAG extraction and multimodal modelling |
| 06 | ONT internal-validation downsampling |
| 07–09 | NGS, EPIC and ULP validation; 07/08 upstream directories precede feature/model steps |
| 10–16 | Cytology comparisons, ROC/CI, genomic landscapes and heatmaps |
| 17 | Performance using pre-existing manual Cluster 2 assignments |
| 18 | Correlations from completed repeated-PoN inputs |
| 19 | Public tissue WGBS full-depth regional methylation/CNA processing |

Module 19 includes a supplied 50-sample manifest. The associated manuscript
reports 136 tissue samples across seven cancer types; the complete 136-sample
manifest and model evaluation were not verified or reconstructed in this release.
The supplied 50-sample material must not be represented as the entire 136-sample
analysis. Tissue labels remain Tumor/Non-malignant control, distinct from PE
mpe/bpe. No WGBS downsampling branch is included. Module 06 is unchanged.

## Input metadata

PE metadata examples retain Sample, group, split and the applicable covariate/
cytology fields. group=mpe/bpe is encoded as 1/0. Original cohort identifiers
must match feature and prediction matrices; EXAMPLE identifiers cannot be joined
to participant data. See metadata/README.md. Do not substitute synthetic values
into manuscript analyses or infer missing patient labels from predictions.

Some scripts reference expected results/ paths that are absent from this release.
These are required external inputs or generated outputs, not bundled data.
Most original Shell templates require output directories to be created first;
module 19 retains its original directory/checkpoint handling.

## Data availability

The associated manuscript reports clinical sequence data in GSA-Human under
HRA019135 and BeadChip IDAT data in OMIX under OMIX017710:

- https://ngdc.cncb.ac.cn/gsa-human/browse/HRA019135
- https://ngdc.cncb.ac.cn/omix/release/OMIX017710

These accessions are transcribed from the manuscript; record visibility, exact
deposited scope and access conditions have not been independently confirmed by
this code-packaging step. Follow the relevant repository's access procedure.
Do not assume that clinical metadata, all derived results or model packages are
part of either deposit. Their release permissions and specific approved access
route remain to be confirmed by the study authors. No access approval or response
deadline is promised by this repository.

The supplied public tissue accessions and annotations are in metadata/public_tissue_wgbs/.
They were inherited from the source package, not independently reverified here.

## Code availability and citation

https://github.com/Che-pathology-lab/nanopore-pe-cfdna-mpe-diagnosis

Use CITATION.cff for the provisional software citation. The final author list,
paper citation and permanent release identifier should be updated when confirmed.
This local package has not been uploaded or assigned an archival DOI.

## License

The analysis code is distributed under the [MIT License](LICENSE), retaining
the existing repository's copyright notice. This license does not grant access
to restricted participant data or override third-party resource licenses.

## Reproducibility limits

The repository includes a synthetic ROC/bootstrap demonstration and its expected
metrics under examples/. Execution scope is described in docs/reproducibility_scope.md.
They do not establish end-to-end clinical reproduction or full journal-policy
compliance. Access to central model artifacts and supporting data for the paper
must be resolved with the authors and the target journal. Known Methods/code
alignment issues are retained in docs/reproducibility_scope.md rather than silently
altering the published analysis. No clinical deployment claim is made by this code release.
