# Input data format

This document specifies the data format expected by `extract_features()`
(`scripts/feature_extractor/feature_extractor.R`), which underlies all three data pathways in this
repository: your own data (`scripts/feature_extractor/run_own_data_template.R`), the synthetic
validation (`scripts/feature_extractor/run_synthetic_validation.R`), and the IBSI phantom
validation (`scripts/feature_extractor/run_ibsi_validation.R`).

## Shape

**Long format: one row per lesion.** Lesions belonging to the same patient
share the same patient identifier and appear as separate rows.

## Required columns

| Column | Argument | Default name | Type | Description |
|---|---|---|---|---|
| Patient identifier | `id_col` | `patient_number` | character or integer | Unique per patient. Repeated across that patient's lesion rows. |
| Uptake value | `value_col` | `SUVmax` | numeric | Lesion-level uptake value. Use `value_col = "SUVmean"` if your data reports SUVmean instead — the extractor doesn't care which, as long as the column name matches what you pass. |

Both are passed as arguments to `extract_features()`, so your data doesn't
have to use these exact column names — just tell the function what you
called them:

```r
extract_features(my_data, category_options = c(4, 6),
                  id_col = "PatientID", value_col = "SUV_mean_bw")
```

**All other columns are ignored.** You do not need to strip your data down
to just these two columns before running the extractor.

## Minimal example

```csv
patient_number,SUVmax
1,3.1
1,4.5
1,2.8
2,6.0
2,5.5
2,7.1
3,2.1
```

Patient 1 has three lesions, patient 2 has three, patient 3 has one. A
worked version of this file, with a few extra (ignored) columns, is at
`data/example/example_lesion_data.csv`.

## The `min_values` parameter

`extract_features(..., min_values = 2)` sets the minimum number of lesions a
patient must have for features to be computed. Patients below that
threshold receive `NA` for every feature (they are **not** dropped from the
output — you still get a row for them, with `number` set to their actual
lesion count and everything else `NA`). Default is 2, since a single lesion
has no meaningful spread, distribution shape, or diversity.

In the minimal example above, with the default `min_values = 2`, patient 3
would get a row of `NA`s.

## Missing values within a patient

**Known limitation:** if a lesion row has `NA` in `value_col`, the extractor
does not currently filter it out automatically, and several of the
underlying statistics (`quantile()` in particular) will error rather than
silently return `NA` for that patient. **Remove or impute rows with missing
uptake values before calling `extract_features()`** — for example:

```r
my_data <- my_data[!is.na(my_data$SUVmax), ]
```

If a patient loses lesions to this filtering and drops below `min_values`,
they will correctly receive an `NA` row rather than erroring.

## Diversity features and `category_options`

`category_options` (e.g. `c(4, 6)`) controls how many uptake categories the
ecological-diversity features (`simpson.D`, `evenness`, `rao.stirling`) are
computed over. Category boundaries are estimated across the **whole cohort**
(via quantile regression), not per patient, so:

- Pass `c()` to skip diversity features entirely — appropriate for
  single-ROI validation (see `run_ibsi_validation.R`, where "diversity
  across lesions" isn't a meaningful concept for a single phantom region).
- A cohort needs enough patients and enough spread in uptake values for the
  category-boundary estimation (`lqmm`) to converge sensibly. There's no
  hard minimum enforced by the code, but very small or very homogeneous
  cohorts may produce unstable category boundaries — inspect the output
  before trusting it if your cohort is small (rule of thumb: tens of
  patients at minimum, not just a handful).

## Output

One row per patient, with:

- `patient_name` — the patient identifier (as a string, regardless of the
  input type)
- `number` — lesion count actually used
- All 22 first-order features (see the feature table in the main README)
- For each value in `category_options`: `simpson.D_<n>`, `evenness_<n>`,
  `rao.stirling_<n>` (e.g. `simpson.D_4`, `evenness_4`, `rao.stirling_4` for
  `category_options = c(4, 6)`)
