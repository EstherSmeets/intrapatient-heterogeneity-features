### Validate intra-patient heterogeneity features against the IBSI digital phantom ###
#
# Validates the first-order features in R/feature_extractor.R against the
# IBSI (Image Biomarker Standardisation Initiative) digital phantom.
#
# Reference: Zwanenburg et al., Radiology 2020; https://doi.org/10.1148/radiol.2020191145
#
# The phantom files (phantom.nii.gz, mask.nii.gz) and the official IBSI
# reference spreadsheet are NOT bundled in this repository -- they are
# Zwanenburg et al.'s materials, not ours to redistribute. Download them
# from the official IBSI repository and place them as described in
# data/ibsi_phantom/README.md before running this script.
#
# Run from the command line:
#   Rscript scripts/run_ibsi_validation.R

#### Libraries ####
library(RNifti)
library(dplyr)
library(tidyr)
library(readxl)
library(openxlsx)
library(here)

source(here::here("scripts/feature_extractor/feature_extractor.R"))

# Variable
img_path   <- here::here("data", "ibsi_phantom", "phantom.nii.gz")
mask_path  <- here::here("data", "ibsi_phantom", "mask.nii.gz")
ref_path   <- here::here("data", "ibsi_phantom", "20181210_lifex.xlsx")
output_dir <- here::here("results")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

category_options <- c() # skip diversity features for phantom validation (single ROI, not multi-lesion)
min_values <- 1 # a single-ROI phantom validation has one "lesion" (the whole ROI) per "patient"

# Import data
img <- readNifti(img_path)
mask <- readNifti(mask_path)
ref_data <- readxl::read_excel(ref_path)

# Select values in mask
roi_vals <- img[mask > 0]
if (length(roi_vals) == 0) stop("No voxels found inside mask. Check IMG_PATH and MASK_PATH.")
cat("IBSI phantom: ", length(roi_vals), "voxels in ROI\n")

# Build single-patient dataset
data_set <- data.frame(
  SUVmax         = roi_vals,
  patient_number = "IBSI_phantom"
)

# Calculate features 
meta_features <- extract_features(data_set, category_options, min_values)

# Import references 
phantom_ref <- ref_data %>%
  filter(family %in% c("Statistics", "Intensity histogram")) %>%
  mutate(image_biomarker = trimws(tolower(image_biomarker)))

name_map <- c(
  energy        = "energy",
  entropy       = "entropy",
  min           = "minimum",
  max           = "maximum",
  mean          = "mean",
  median        = "median",
  variance      = "variance",
  skewness      = "skewness",
  kurtosis      = "kurtosis",
  range         = "range",
  iqr           = "interquartile range",
  mad           = "mean absolute deviation",
  robust_mad    = "robust mean absolute deviation",
  rms           = "root mean square",
  cov           = "coefficient of variation"
)

# IBSI match logic
ibsi_match <- function(computed, reference, tolerance) {
  if (is.na(computed) || is.na(reference)) return(NA_character_)
  if (reference == 0) {
    diff_val <- abs(computed - reference)
    tol_eff  <- if (tolerance > 0) tolerance else 0.0005
  } else {
    n_dec    <- 3 - 1 - floor(log10(abs(reference)))
    diff_val <- abs(round(computed - reference, n_dec))
    tol_eff  <- if (tolerance > 0) tolerance else 10^(-n_dec) / 2
  }
  if (diff_val <= tol_eff)           "match"
  else if (diff_val < 3 * tol_eff)  "partial match"
  else                               "no match"
}

# Builed comparison table 
my_long <- meta_features %>%
  pivot_longer(-patient_name, names_to = "feature_code", values_to = "computed_value") %>%
  mutate(ibsi_name = name_map[feature_code]) %>%
  filter(!is.na(ibsi_name)) %>%
  mutate(ibsi_name = trimws(tolower(ibsi_name)))

# Primary: Statistics family. For entropy: use Intensity histogram reference.
ref_stats <- phantom_ref %>%
  filter(family == "Statistics") %>%
  select(image_biomarker, family, benchmark_value, tolerance, tag)

ref_ih_entropy <- phantom_ref %>%
  filter(family == "Intensity histogram", image_biomarker == "entropy") %>%
  select(image_biomarker, family, benchmark_value, tolerance, tag)

ref_combined <- bind_rows(ref_stats, ref_ih_entropy) %>%
  distinct(image_biomarker, .keep_all = TRUE)

compare_df <- my_long %>%
  left_join(ref_combined, by = c("ibsi_name" = "image_biomarker")) %>%
  rowwise() %>%
  mutate(
    match_status = ibsi_match(computed_value, benchmark_value, tolerance),
    difference   = round(abs(computed_value - benchmark_value), 6)
  ) %>%
  ungroup() %>%
  arrange(family, ibsi_name) %>%
  select(family, ibsi_name, benchmark_value, tolerance,
         computed_value, difference, match_status, tag)

cat("\n--- IBSI comparison ---\n")
print(as.data.frame(compare_df), digits = 6)

n_match   <- sum(compare_df$match_status == "match",         na.rm = TRUE)
n_partial <- sum(compare_df$match_status == "partial match", na.rm = TRUE)
n_no      <- sum(compare_df$match_status == "no match",      na.rm = TRUE)
n_total   <- sum(!is.na(compare_df$match_status))
cat(sprintf(
  "\nResult: %d match | %d partial match | %d no match  (of %d features with IBSI reference)\n",
  n_match, n_partial, n_no, n_total))

# Write computed values into a copy of the IBSI xlsx
wb         <- loadWorkbook(ref_path)
sheet_name <- "digital phantom"

ref_raw  <- read_excel(ref_path, sheet = sheet_name)
res_col  <- which(colnames(ref_raw) == "your_result")

for (i in seq_len(nrow(compare_df))) {
  row_tag <- compare_df$tag[i]
  if (is.na(row_tag)) next
  xlsx_row <- which(ref_raw$tag == row_tag)
  if (length(xlsx_row) == 0) next
  writeData(wb,
            sheet    = sheet_name,
            x        = compare_df$computed_value[i],
            startCol = res_col,
            startRow = xlsx_row + 1,   # +1 for header row
            colNames = FALSE)
}

out_xlsx <- file.path(output_dir, "IBSI_phantom_results.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)


