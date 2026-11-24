# intrapatient-heterogeneity-features
Extract intrapatient heterogeneity features from PET-features of metastatic patients

R implementation of a framework for quantifying intra-patient heterogeneity from lesion-level PET uptake values (SUVmax or SUVmean). 
Computes 22 features across eight conceptual groups: central tendency, spread, distributional shape, entropy/energy, statistical inequality, and distance-weighted ecological diversity — from any dataset of per-lesion uptake values.

This repository accompanies the methods paper *[Paper title ]* (... et al., in preparation), and provides:

- the feature extraction function itself, for use on your own lesion-level PET data;
- the synthetic data generator used for the paper's systematic validation (fully reproducible, seed = 42);
- validation against the [IBSI digital phantom](https://doi.org/10.1148/radiol.2020191145) for the first-order feature subset.

## Contents

- [Feature groups](#feature-groups)
- [Installation](#installation)
- [Repository structure](#repository-structure)
- [Usage](#usage)
  - [A. Extract features from your own data](#a-extract-features-from-your-own-data)
  - [B. Reproduce the synthetic validation](#b-reproduce-the-synthetic-validation)
  - [C. Validate against the IBSI digital phantom](#c-validate-against-the-ibsi-digital-phantom)
- [Input data format](#input-data-format)
- [Output](#output)
- [Reproducibility](#reproducibility)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## Feature groups

| Group | Features |
|---|---|
| Central tendency | `mean`, `median`, `rms`, `energy_normalized` |
| Minimum | `min` |
| Spread | `max`, `variance`, `ipr`, `mad`, `iqr`, `robust_mad` |
| Distributional shape | `skewness`, `kurtosis` |
| Entropy and energy | `entropy`, `energy` |
| Statistical inequality | `gini`, `theil`, `cov`, `atkinson_low`, `atkinson_high`, `hoover` |
| Ecological diversity | `simpson.D_4`, `evenness_4`, `simpson.D_6`, `evenness_6` |
| Distance-weighted diversity | `rao.stirling_4`, `rao.stirling_6` |

Feature redundancy and a recommended representative subset (one feature per group, selected by hierarchical clustering on absolute Spearman correlation) are described in the accompanying paper. 
`gini` is recommended as the representative for the inequality group; `min` should always be retained as an independent feature regardless of clustering outcome.

## Installation

Requires R ≥ 4.x. Install dependencies:

```r
install.packages(c("dplyr", "tidyr", "ineq", "lqmm", "here",
                    "MASS", "readxl", "openxlsx", "RNifti"))
```

Clone the repository:

```bash
git clone https://github.com/[YOUR-GITHUB-USERNAME]/[REPO-NAME].git
cd [REPO-NAME]
```

## Repository structure

```
├── scripts/
│   ├── feature_extractor/
        ├── feature_extractor.R        # core extract_features() function
│       ├── run_synthetic_validation.R # synthetic data generation + feature extraction
│       ├── run_own_data_template.R    # template for applying the extractor to your own data
│       └── run_ibsi_validation.R      # validation against the IBSI digital phantom
├── data/
│   ├── example/                   # small fabricated example CSVs (schema demonstration only)
│   └── ibsi_phantom/              # download instructions for the official IBSI phantom files
├── docs/
│   └── data_format.md             # full input data format specification
└── results/                       # script output (created on run, not version-controlled)
```

## Usage

All scripts are run from the command line via `Rscript` and use `here::here()` for paths, so run them from anywhere inside the cloned repository.

### A. Extract features from your own data

`scripts/feature_extractor/run_own_data_template.R` is a template, not a plug-and-play script — copy it and adapt the column-selection lines to your own dataset's column names, then point it at your CSV. 
See [Input data format](#input-data-format) below for the schema the extractor itself expects after your adaptation.

```bash
Rscript scripts/feature_extractor/run_own_data_template.R --input path/to/your_lesion_data.csv --output results/
```

A minimal fabricated example is provided at `data/example/example_lesion_data.csv` so you can confirm the pipeline runs before pointing it at real data.

### B. Reproduce the synthetic validation

Regenerates the 1×3×5 fully factorial synthetic validation described in the paper (1 distributions × 3 uptake range scenarios × 5 lesion-count groups, 1,000 synthetic patients per condition). Fixed at `seed = 42`, so output is deterministic.
Both normal and lognormal distributions are an option.

```bash
Rscript scripts/feature_extractor/run_syntetic_validation.R
```

Output: `results/Meta_features_<n_lesions>.csv` and `results/PET_data_cohort_<n_lesions>.csv` per lesion-count group.

### C. Validate against the IBSI digital phantom

Reproduces the first-order feature validation against the [IBSI digital phantom](https://doi.org/10.1148/radiol.2020191145) reference values. 
Requires the official phantom files, which are **not redistributed in this repository** — see `data/ibsi_phantom/README.md` for the download link and expected file placement.

# IBSI digital phantom

`scripts/feature_extractor/run_ibsi_validation.R` validates the first-order features against the
[IBSI digital phantom](https://doi.org/10.1148/radiol.2020191145) (Zwanenburg et al., Radiology 2020).

These files are **not included in this repository** — they belong to the IBSI
project, not to us. Download them from the official IBSI repository and place
them here as:

```
data/ibsi_phantom/
├── phantom.nii.gz
├── mask.nii.gz
└── 20181210_lifex.xlsx
```

```bash
Rscript scripts/feature_extractor/run_ibsi_validation.R
```

## Input data format

The extractor (`extract_features()`) expects a **long-format data frame**: one row per lesion, with lesions belonging to the same patient sharing the same identifier.

| Column | Required | Type | Description |
|---|---|---|---|
| `patient_number` | Yes | character or integer | Unique patient identifier. |
| `SUVmax` | Yes | numeric | Lesion-level uptake value. *(Note: the current implementation reads this column under the literal name `SUVmax`; if your data uses `SUVmean`, rename the column before calling the extractor — see `docs/data_format.md` for details.)* |

Any additional columns in your data (tracer, lesion site, scan date, etc.) are ignored by `extract_features()` and do not need to be removed.

Patients with fewer lesions than `min_values` (a function argument) receive `NA` for all features rather than being silently dropped.

Full schema documentation, including the diversity-feature category options and worked examples, is in [`docs/data_format.md`](docs/data_format.md).

## Output

`extract_features()` returns one row per patient with:

- `patient_name`, `number` (lesion count used)
- all first-order features listed in [Feature groups](#feature-groups)
- diversity features for each requested category count (e.g. `simpson.D_4`, `evenness_4`, `rao.stirling_4` for `category_options = 4`)

## Reproducibility

- All random processes use `set.seed(42)`.
- Paths are resolved with `here::here()` — no hardcoded absolute paths, no RStudio-only path resolution.
- Package versions are pinned via `renv.lock` (see Installation).

## Citation

If you use this code, please cite:

> Smeets E.M.M., Noortman W.A.  *[Paper title]*. Manuscript in preparation, [target journal].

A `CITATION.cff` file is provided; GitHub will render a "Cite this repository" option automatically. 
This entry will be updated with the full citation and DOI once the paper is published, and a Zenodo-archived, versioned release of the code will be linked here at that point.

## License

[MIT](LICENSE)

## Contact

Esther Smeets, PostDoc, University Medical Center Utrecht
Wyanne Noortman, PostDoc, UMCG Groningen
Supervisors: Sjoerd Elias

For questions or issues, please open a [GitHub issue](https://github.com/[YOUR-GITHUB-USERNAME]/[REPO-NAME]/issues) rather than emailing directly, so answers are visible to other users.
