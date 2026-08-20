#### Calculate intra-patient heterogeneity features from synthetic data ####
#
# Generates synthetic lesion-level uptake data and computes intra-patient
# heterogeneity features across a grid of uptake-range scenarios and lesion
# counts, for ONE distribution family per run.
#
# Per the project's validation design, SUVmax uses a log-normal distribution
# and SUVmean uses a normal distribution -- these are never mixed in a single
# run. Run this script twice (once with `distribution <- "lognormal"`, once
# with `distribution <- "normal"`) to cover both.
#
# Run from the command line:
#   Rscript scripts/run_synthetic_validation.R

library(MASS)
library(here)

source(here::here("scripts/feature_extractor/feature_extractor.R"))

set.seed(42)

#### Configuration ####

# Distribution family for this run.
#   "lognormal" -> values generated on the log scale, exp()-transformed back
#                  (use for SUVmax)
#   "normal"    -> values generated directly on the linear scale
#                  (use for SUVmean)
distribution <- "lognormal"
value_col    <- if (distribution == "lognormal") "SUVmax" else "SUVmean"

n_patients <- 10
lesions_per_patient_list <- c(5, 10, 20, 40, 80)
category_options <- c(4, 6)
min_values <- 2

# Uptake-range scenarios, one row per clinical scenario cited in the paper.
# `clinical_range` is the literature-reported range in LINEAR SUVmax units --
# this is the only part of each scenario anchored to a citation.
#
# `mean_bounds` are derived from `clinical_range` below (log() for lognormal,
# a direct linear scaling for normal/SUVmean), NOT specified directly, so
# that a distribution's sampled values actually land near the cited clinical
# range once transformed back to linear units. Getting this derivation wrong
# (i.e. treating clinical_range as if it were already log-scale, or applying
# exp() without first taking the log going in) is what produced physically
# implausible values (SUVmax in the hundreds) when this was tested.
#
# `sd_within_patient` (log-scale SD if lognormal, linear SD if normal) and
# `correlation_bounds` are simulation design choices, not literature values.
# The defaults below are a reasonable starting point (chosen to keep the 5th-
# 95th percentile of simulated values roughly within clinical_range) but you
# should check the printed diagnostic quantiles against your expectations
# before treating results as final -- these govern how *spread out* lesion
# uptake is within a simulated patient, which is exactly what the paper's
# heterogeneity features are measuring, so getting this right matters more
# than almost any other parameter here.
range_scenarios <- list(
  low = list(
    # Breast HR+; Aras et al. 2008, Eur J Radiol: SUVmax range 1.5-5
    clinical_range      = c(1.5, 5),
    sd_within_patient    = c(0.10, 0.35),
    correlation_bounds  = c(0, 0.9)
  ),
  intermediate = list(
    # Metastatic solid tumours; Vogsen et al. 2018 (PMC5943993): SUVmax range 2-10
    clinical_range      = c(2, 10),
    sd_within_patient    = c(0.10, 0.35),
    correlation_bounds  = c(0, 0.9)
  ),
  high = list(
    # R/R DLBCL; Nols et al. 2011 (PMID 21512731) / Persky et al. 2010 ASCO: SUVmax range 2-25
    clinical_range      = c(2, 25),
    sd_within_patient    = c(0.10, 0.35),
    correlation_bounds  = c(0, 0.9)
  )
)

# Nahmias & Wahl 2008, JNM 49(11): SUVmean/SUVmax ratio ~ 0.52-0.55, used to
# scale clinical_range down when distribution == "normal" (SUVmean).
suvmean_suvmax_ratio <- c(0.52, 0.55)

#' Derive mean_bounds for a scenario, on the scale the sampling distribution expects
#'
#' @param clinical_range Length-2 linear-scale SUVmax range from the literature.
#' @param distribution "lognormal" or "normal".
#' @return Length-2 numeric vector: log-scale bounds if lognormal, linear
#'   SUVmean-scaled bounds if normal.
derive_mean_bounds <- function(clinical_range, distribution) {
  if (distribution == "lognormal") {
    log(clinical_range)
  } else {
    clinical_range * mean(suvmean_suvmax_ratio)
  }
}

output_dir <- here::here("results")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#' Simulate one lesion-level cohort for a given scenario
#'
#' Draws per-patient mean, SD, and intra-patient correlation from uniform
#' distributions bounded by \code{mean_bounds}/\code{sd_bounds}/
#' \code{correlation_bounds}, then simulates \code{n_lesions} correlated
#' lesion-level values per patient from a multivariate normal distribution.
#' If \code{distribution == "lognormal"}, the sampled values are treated as
#' log-scale and exponentiated before being returned.
#'
#' @param n_patients Number of synthetic patients to simulate.
#' @param n_lesions Number of lesions per patient (fixed across patients in
#'   a single call).
#' @param mean_bounds Length-2 numeric vector, uniform sampling bounds for
#'   each patient's mean (log-scale mean if \code{distribution == "lognormal"}).
#' @param sd_bounds Length-2 numeric vector, uniform sampling bounds for each
#'   patient's SD (log-scale SD if \code{distribution == "lognormal"}).
#' @param correlation_bounds Length-2 numeric vector, uniform sampling bounds
#'   for each patient's intra-patient lesion correlation.
#' @param distribution Either \code{"lognormal"} or \code{"normal"}.
#' @param scenario_label Character label stored in the output for traceability
#'   (e.g. \code{"low"}, \code{"intermediate"}, \code{"high"}).
#' @param value_col Name to give the simulated-value column in the output
#'   (\code{"SUVmax"} or \code{"SUVmean"}).
#'
#' @return A long-format data frame with columns \code{patient_number},
#'   \code{value_col}, \code{correlation}, and \code{scenario}.
simulate_lesion_cohort <- function(n_patients, n_lesions,
                                   mean_bounds, sd_bounds, correlation_bounds,
                                   distribution, scenario_label, value_col) {
  
  correlations <- stats::runif(n_patients, correlation_bounds[1], correlation_bounds[2])
  means        <- stats::runif(n_patients, mean_bounds[1], mean_bounds[2])
  sds          <- stats::runif(n_patients, sd_bounds[1], sd_bounds[2])
  
  patient_rows <- vector("list", n_patients)
  
  for (i in seq_len(n_patients)) {
    corr_matrix <- matrix(correlations[i], nrow = n_lesions, ncol = n_lesions)
    diag(corr_matrix) <- 1
    cov_matrix <- (sds[i]^2) * corr_matrix
    mean_vec   <- rep(means[i], n_lesions)
    
    sampled <- t(MASS::mvrnorm(1, mu = mean_vec, Sigma = cov_matrix))
    
    values <- if (distribution == "lognormal") exp(sampled) else sampled
    
    patient_id <- paste0(scenario_label, "_", n_lesions, "_", i)
    patient_rows[[i]] <- data.frame(
      patient_number = patient_id,
      value          = as.numeric(values),
      correlation    = correlations[i],
      scenario       = scenario_label
    )
  }
  
  out <- do.call(rbind, patient_rows)
  names(out)[names(out) == "value"] <- value_col
  out
}

#### Run simulation across all range scenarios and lesion counts ####

for (scenario_name in names(range_scenarios)) {
  scenario <- range_scenarios[[scenario_name]]
  cat(sprintf("\n=== Scenario: %s (%s) ===\n", scenario_name, distribution))
  
  mean_bounds <- derive_mean_bounds(scenario$clinical_range, distribution)
  
  for (n_lesions in lesions_per_patient_list) {
    cat("  Simulating", n_lesions, "lesions per patient\n")
    
    data_set <- simulate_lesion_cohort(
      n_patients         = n_patients,
      n_lesions           = n_lesions,
      mean_bounds         = mean_bounds,
      sd_bounds            = scenario$sd_within_patient,
      correlation_bounds  = scenario$correlation_bounds,
      distribution         = distribution,
      scenario_label       = scenario_name,
      value_col            = value_col
    )
    
    # Outlier trimming (99.99th percentile cap, floor at 0)
    p9999 <- stats::quantile(data_set[[value_col]], 0.9999, na.rm = TRUE)
    data_set[[value_col]] <- pmin(data_set[[value_col]], p9999)
    data_set[[value_col]] <- pmax(data_set[[value_col]], 0)
    
    # Diagnostic: verify simulated values actually land near clinical_range.
    # Check this against the cited literature range before trusting results.
    q <- stats::quantile(data_set[[value_col]], c(.05, .5, .95))
    cat(sprintf("    Achieved %s range (P5-P50-P95): %.2f - %.2f - %.2f  |  target clinical_range: %.1f-%.1f\n",
                value_col, q[1], q[2], q[3], scenario$clinical_range[1], scenario$clinical_range[2]))
    
    meta_features <- extract_features(
      data_set, category_options, min_values,
      value_col = value_col
    )
    
    write.csv(
      meta_features,
      file.path(output_dir, sprintf("Meta_features_%s_%s_%d.csv", distribution, scenario_name, n_lesions)),
      row.names = FALSE
    )
    write.csv(
      data_set,
      file.path(output_dir, sprintf("PET_data_cohort_%s_%s_%d.csv", distribution, scenario_name, n_lesions)),
      row.names = FALSE
    )
  }
}

cat("\nDone. Output written to", output_dir, "\n")