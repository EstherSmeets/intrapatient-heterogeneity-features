#### Intra-patient heterogeneity feature extractor ####
#
# Core feature computation used by all three data pathways:
#   - scripts/run_synthetic_validation.R
#   - scripts/run_own_data_template.R
#   - scripts/run_ibsi_validation.R
#
# This file defines functions only; it is meant to be sourced, not run directly.

library(ineq)
library(lqmm)

#' Shannon entropy via epsilon-corrected fixed binning over a global range
#'
#' Alternative entropy estimator that bins values into 10 equal-width bins
#' spanning a caller-supplied global range (e.g. the full cohort's min/max),
#' rather than a per-patient range. Not currently used by \code{extract_features}
#' (which uses \code{compute_entropy_fixed_width} instead) but kept for
#' reference / comparison, since binning against a shared global range can be
#' preferable when comparing entropy across patients with very different value
#' ranges.
#'
#' @param values Numeric vector of lesion-level uptake values for one patient.
#' @param global_min_suv Minimum value across the full cohort, used as the
#'   lower bin edge.
#' @param global_max_suv Maximum value across the full cohort, used as the
#'   upper bin edge.
#' @param epsilon Small constant added inside the log to avoid log(0) for
#'   empty bins.
#' @return Numeric scalar: Shannon entropy in bits.
compute_entropy_epsilon <- function(values, global_min_suv, global_max_suv, epsilon = 1e-12) {
  # NOTE: bin edges must run low -> high for cut() to produce valid bins.
  breaks <- seq(global_min_suv, global_max_suv, length.out = 11)
  bins <- cut(values, breaks = breaks, include.lowest = TRUE, right = FALSE)
  
  counts <- table(bins)
  probs <- counts / sum(counts)
  
  entropy <- -sum(probs * log2(probs + epsilon), na.rm = TRUE)
  return(entropy)
}

#' Shannon entropy via fixed-width binning
#'
#' Bins a patient's values into fixed-width bins anchored at \code{anchor},
#' spanning from below the patient's minimum to above the patient's maximum,
#' and computes Shannon entropy over the resulting bin-occupancy distribution.
#' This is the entropy estimator used by \code{extract_features}.
#'
#' @param values Numeric vector of lesion-level uptake values for one patient.
#' @param binwidth Width of each bin, in the same units as \code{values}
#'   (e.g. SUV units). Passed through from \code{extract_features}'s
#'   \code{entropy_binwidth} argument.
#' @param anchor Bin edge alignment point. Default 0.
#' @return Numeric scalar: Shannon entropy in bits, or \code{NA_real_} if
#'   fewer than 2 non-missing values are supplied.
compute_entropy_fixed_width <- function(values, binwidth, anchor = 0) {
  v <- values[!is.na(values)]
  if (length(v) <= 1 || binwidth <= 0) return(NA_real_)
  
  vmax <- max(v)
  vmin <- min(v)
  
  breaks <- seq(anchor, vmax + binwidth, by = binwidth)
  if (vmin < anchor) {
    extra <- seq(anchor - binwidth,
                 floor((vmin - anchor) / binwidth) * binwidth + anchor,
                 by = -binwidth)
    breaks <- sort(unique(c(breaks, extra)))
  }
  
  bins <- cut(v, breaks = breaks, include.lowest = TRUE, right = FALSE)
  counts <- as.numeric(table(bins))
  probs <- counts[counts > 0] / sum(counts)
  
  -sum(probs * log2(probs))
}

#' Ecological diversity metrics for one patient's categorical distribution
#'
#' Computes Simpson's diversity index, evenness, and a distance-weighted
#' Rao-Stirling diversity index from a categorical probability distribution
#' and a between-category distance matrix.
#'
#' @param p Numeric vector of category proportions (need not already sum to
#'   1; will be renormalised). Negative values are clipped to 0.
#' @param distance_matrix A k x k numeric matrix of between-category
#'   distances, where k == length(p).
#' @return A list with elements \code{simpson.D}, \code{evenness}, and
#'   \code{rao.stirling}.
compute_diversity_metrics <- function(p, distance_matrix) {
  p <- as.numeric(p)
  p[p < 0] <- 0
  s <- sum(p)
  if (s > 0) p <- p / s
  
  k <- length(p)
  
  simpson_D <- if (k > 1) 1 - sum(p^2) else 0
  evenness  <- if (k > 1) simpson_D / (1 - 1 / k) else NA_real_
  
  if (!is.matrix(distance_matrix) ||
      nrow(distance_matrix) != k ||
      ncol(distance_matrix) != k) {
    stop("distance_matrix must be a k x k matrix where k == length(p)")
  }
  
  rao <- as.numeric(t(p) %*% distance_matrix %*% p)
  
  list(
    simpson.D    = simpson_D,
    evenness     = evenness,
    rao.stirling = rao
  )
}

#' Extract intra-patient heterogeneity features from lesion-level uptake data
#'
#' Computes 22 first-order and distribution-shape features (central
#' tendency, spread, distributional shape, entropy/energy, and statistical
#' inequality) plus ecological-diversity features (Simpson diversity,
#' evenness, Rao-Stirling diversity) for each patient in \code{data_set},
#' from that patient's set of lesion-level uptake values.
#'
#' @param data_set A data frame in long format: one row per lesion. Must
#'   contain a patient identifier column (\code{id_col}) and a numeric
#'   uptake-value column (\code{value_col}). Extra columns are ignored.
#' @param category_options Integer vector of category counts to compute
#'   diversity features for (e.g. \code{c(4, 6)}). Pass \code{c()} to skip
#'   diversity features entirely (e.g. for single-ROI phantom validation,
#'   where diversity across "lesions" is not meaningful).
#' @param min_values Minimum number of lesions a patient must have for
#'   features to be computed; patients with fewer lesions receive \code{NA}
#'   for all features. Default 2 (a single lesion has no meaningful spread).
#' @param value_col Name of the column in \code{data_set} holding the
#'   uptake value to use (e.g. \code{"SUVmax"} or \code{"SUVmean"}).
#'   Default \code{"SUVmax"}.
#' @param id_col Name of the column in \code{data_set} holding the patient
#'   identifier. Default \code{"patient_number"}.
#' @param entropy_binwidth Bin width (in \code{value_col} units) used by
#'   \code{compute_entropy_fixed_width}. Default 0.5.
#'
#' @return A data frame with one row per patient, containing \code{patient_name},
#'   \code{number} (lesion count used), the first-order features, and — for
#'   each entry in \code{category_options} — \code{simpson.D_<n>},
#'   \code{evenness_<n>}, \code{rao.stirling_<n>}.
#'
#' @examples
#' \dontrun{
#'   df <- data.frame(
#'     patient_number = c(1, 1, 1, 2, 2),
#'     SUVmax         = c(3.1, 4.5, 2.8, 6.0, 5.5)
#'   )
#'   extract_features(df, category_options = c(4), min_values = 2)
#' }
extract_features <- function(data_set,
                             category_options,
                             min_values = 2,
                             value_col = "SUVmax",
                             id_col = "patient_number",
                             entropy_binwidth = 0.5) {
  
  if (!value_col %in% names(data_set)) {
    stop(sprintf("value_col '%s' not found in data_set. Available columns: %s",
                 value_col, paste(names(data_set), collapse = ", ")))
  }
  if (!id_col %in% names(data_set)) {
    stop(sprintf("id_col '%s' not found in data_set. Available columns: %s",
                 id_col, paste(names(data_set), collapse = ", ")))
  }
  
  values_all <- data_set[[value_col]]
  ids_all    <- data_set[[id_col]]
  
  # Global range, used only by compute_entropy_epsilon (not the default path).
  global_min_suv <- min(values_all, na.rm = TRUE)
  global_max_suv <- max(values_all, na.rm = TRUE)
  
  # Canonical, ordered list of first-order feature names. Both the "enough
  # lesions" branch and the "too few lesions" (NA) branch below MUST produce
  # a named vector with exactly these names, in this order — this list is
  # the single source of truth that keeps the two branches in sync.
  feature_names <- c(
    "number", "energy", "energy_normalized", "entropy",
    "min", "max", "mean", "median", "variance",
    "skewness", "kurtosis", "iqr", "ipr", "mad", "robust_mad",
    "rms", "cov", "gini", "theil", "atkinson_low", "atkinson_high", "hoover"
  )
  
  patient_list <- unique(ids_all)
  meta_features <- vector("list", length(patient_list))
  names(meta_features) <- as.character(patient_list)
  
  for (patient in patient_list) {
    
    values <- values_all[ids_all == patient]
    n <- length(values)
    
    if (n < min_values) {
      row <- stats::setNames(rep(NA_real_, length(feature_names)), feature_names)
      row["number"] <- n
      meta_features[[as.character(patient)]] <- row
      next
    }
    
    max_suv <- max(values)
    min_suv <- min(values)
    mu      <- mean(values)
    dev     <- values - mu
    m2      <- mean(dev^2)
    sdpop   <- sqrt(m2)
    
    qs  <- stats::quantile(values, c(.10, .25, .50, .75, .90), names = FALSE)
    p10 <- qs[1]; q1 <- qs[2]; med <- qs[3]; q3 <- qs[4]; p90 <- qs[5]
    
    values_1090 <- values[values >= p10 & values <= p90]
    mu_1090 <- if (length(values_1090) > 0) mean(values_1090) else NA_real_
    
    all_pos <- all(values > 0)
    energy  <- sum(values^2)
    
    row <- c(
      number            = n,
      energy            = energy,
      energy_normalized = energy / n,
      entropy           = compute_entropy_fixed_width(values, entropy_binwidth),
      min               = min_suv,
      max               = max_suv,
      mean              = mu,
      median            = med,
      variance          = m2,
      skewness          = if (sdpop > 0) mean((dev / sdpop)^3) else 0,
      kurtosis          = if (sdpop > 0) mean((dev / sdpop)^4) - 3 else 0,
      iqr               = q3 - q1,
      ipr               = p90 - p10,
      mad               = mean(abs(values - mu)),
      robust_mad        = if (length(values_1090) > 0) mean(abs(values_1090 - mu_1090)) else NA_real_,
      rms               = sqrt(mean(values^2)),
      cov               = if (mu != 0) sdpop / mu else NA_real_,
      # Gini is defined for non-negative values (uptake values are always
      # >= 0), so it is computed unconditionally, unlike Theil/Atkinson/CoV
      # which require strictly positive values (they involve logs/division).
      gini              = ineq::Gini(values),
      theil             = if (all_pos) ineq::Theil(values) else NA_real_,
      atkinson_low      = if (all_pos) ineq::Atkinson(values, parameter = 0.4) else NA_real_,
      atkinson_high     = if (all_pos) ineq::Atkinson(values, parameter = 1.8) else NA_real_,
      hoover            = if (all_pos) ineq::RS(values) else NA_real_
    )
    
    # Defensive check: catch any future edit that adds/removes/reorders a
    # feature in one branch without updating the other, before it can
    # silently corrupt output the way the original bug did.
    stopifnot(identical(names(row), feature_names))
    
    meta_features[[as.character(patient)]] <- row
  }
  
  # IMPORTANT: bind_rows() (not do.call(rbind, ...)) aligns columns by name.
  # do.call(rbind, ...) on named vectors of differing length/order silently
  # mis-aligns values instead of raising an error — this was the source of
  # a critical data-corruption bug in an earlier version of this function.
  feature_df <- dplyr::bind_rows(lapply(meta_features, as.list))
  feature_df$patient_name <- names(meta_features)
  
  # Compute diversity features
  selection_combined <- data.frame(
    patient_number = as.character(unique(ids_all))
  )
  names(selection_combined) <- id_col
  
  for (n_cat in category_options) {
    
    cat("Number of categories:", n_cat, "\n")
    
    tau_vec <- seq(0, 1, length.out = n_cat + 1)[2:n_cat]  # interior quantiles
    
    # lqmm stores the unevaluated call and re-evaluates parts of it later
    # (e.g. when summary() builds the boundary table), in an environment
    # where local variables from this function (like `id_col`) no longer
    # exist. Passing group = data_set[[id_col]] directly works at fit time
    # but breaks on that later re-evaluation ("object 'id_col' not found").
    # Fix: give lqmm a literal column name to reference instead, by adding
    # it to a local copy of the data under a fixed name.
    fit_data <- data_set
    fit_data$.lqmm_group <- fit_data[[id_col]]
    
    set.seed(42)
    fit_formula <- stats::as.formula(paste(value_col, "~ 1"))
    fit <- lqmm::lqmm(
      fixed = fit_formula, random = ~1, group = .lqmm_group,
      data = fit_data, tau = tau_vec, nK = 11, type = "normal",
      control = lqmm::lqmmControl(method = "df", UP_max_iter = 500)
    )
    
    summary_table <- summary(fit)$tTable
    boundary_values <- vapply(
      as.character(tau_vec),
      function(tau) summary_table[[tau]][, 1],
      numeric(1)
    )
    breaks_cat <- c(-Inf, boundary_values, Inf)
    
    category <- as.integer(
      cut(values_all, breaks = breaks_cat, labels = FALSE, include.lowest = TRUE)
    )
    
    # Distance matrix: SUV values with larger differences get larger distance.
    finite_breaks <- boundary_values
    midpoints <- c(
      mean(c(global_min_suv, if (length(finite_breaks) > 0) finite_breaks[1] else global_max_suv)),
      if (length(finite_breaks) > 1)
        (finite_breaks[-length(finite_breaks)] + finite_breaks[-1]) / 2,
      mean(c(if (length(finite_breaks) > 0) finite_breaks[length(finite_breaks)] else global_min_suv,
             global_max_suv))
    )
    distance_matrix <- as.matrix(stats::dist(midpoints, method = "euclidean"))
    
    div_rows <- lapply(patient_list, function(pt) {
      cats_pt <- category[ids_all == pt]
      cats_pt <- cats_pt[!is.na(cats_pt)]
      
      if (length(cats_pt) == 0) {
        return(data.frame(
          patient_number = as.character(pt),
          simpson.D      = NA_real_,
          evenness       = NA_real_,
          rao.stirling   = NA_real_
        ))
      }
      
      counts <- as.numeric(table(factor(cats_pt, levels = 1:n_cat)))
      total  <- sum(counts)
      
      if (total == 0) {
        return(data.frame(
          patient_number = as.character(pt),
          simpson.D      = NA_real_,
          evenness       = NA_real_,
          rao.stirling   = NA_real_
        ))
      }
      
      p_vec <- counts / total
      metrics <- compute_diversity_metrics(p_vec, distance_matrix)
      
      data.frame(
        patient_number = as.character(pt),
        simpson.D      = metrics$simpson.D,
        evenness       = metrics$evenness,
        rao.stirling   = metrics$rao.stirling
      )
    })
    
    diversity_out <- dplyr::bind_rows(div_rows)
    names(diversity_out)[1] <- id_col
    colnames(diversity_out)[2:4] <- paste0(colnames(diversity_out)[2:4], "_", n_cat)
    
    selection_combined <- dplyr::left_join(selection_combined, diversity_out, by = id_col)
  }
  
  # Join first-order features with diversity features
  feature_df <- dplyr::left_join(
    feature_df, selection_combined,
    by = stats::setNames(id_col, "patient_name")
  )
  
  # Reorder: patient_name and number first
  first_col  <- "patient_name"
  second_col <- "number"
  other_cols <- setdiff(colnames(feature_df), c(first_col, second_col))
  feature_df <- feature_df[, c(first_col, second_col, other_cols)]
  
  return(feature_df)
}