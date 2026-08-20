#### Extract intra-patient heterogeneity features from your own data ####
#
# TEMPLATE SCRIPT -- this is not plug-and-play. Copy it, then adapt the
# "ADAPT THIS" block below to your dataset's actual column names before
# running. See docs/data_format.md for the full input schema.
#
# Run from the command line, e.g.:
#   Rscript scripts/run_own_data_template.R --input data/my_lesion_data.csv --output results/
#   Rscript scripts/run_own_data_template.R --input data/my_lesion_data.csv --value-col SUVmean

library(here)

source(here::here("scripts/feature_extractor/feature_extractor.R"))

#### Parse command-line arguments ####

#' Parse simple --flag value command-line arguments
#'
#' @param args Character vector, typically from commandArgs(trailingOnly = TRUE).
#' @return A named list of argument values.
parse_args <- function(args) {
  defaults <- list(
    input      = NULL,
    output     = here::here("results"),
    `value-col` = "SUVmax",
    `min-values` = "2"
  )
  i <- 1
  while (i <= length(args)) {
    flag <- sub("^--", "", args[i])
    if (flag %in% names(defaults)) {
      defaults[[flag]] <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  defaults
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (is.null(args$input)) {
  stop("Usage: Rscript run_own_data_template.R --input path/to/data.csv [--output results/] [--value-col SUVmax] [--min-values 2]")
}

value_col  <- args[["value-col"]]
min_values <- as.integer(args[["min-values"]])
output_dir <- args$output
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

category_options <- c(4, 6)

#### Import data ####

data_set <- read.csv(args$input)

#### ADAPT THIS: map your columns to the required schema ####
#
# extract_features() requires two columns (see docs/data_format.md):
#   - patient_number : unique patient identifier
#   - <value_col>     : lesion-level uptake value (SUVmax or SUVmean)
#
# Example adaptation, shown commented out, for a clinical export with
# different column names and a tracer filter -- replace with your own:
#
# data_set <- data_set[data_set$tracer == "FDG", ]
# data_set <- data_set[!is.na(data_set$tracer), ]
# names(data_set)[names(data_set) == "PatientID"] <- "patient_number"
# names(data_set)[names(data_set) == "SUV_max_bw"] <- value_col

if (!"patient_number" %in% names(data_set) || !value_col %in% names(data_set)) {
  stop(sprintf(
    paste("data_set must contain 'patient_number' and '%s' columns after adaptation.",
          "Found columns: %s. Edit the 'ADAPT THIS' block above to map your own",
          "column names to these two required columns."),
    value_col, paste(names(data_set), collapse = ", ")
  ))
}

#### Extract intra-patient heterogeneity features ####

meta_features <- extract_features(data_set, category_options, min_values, value_col = value_col)

#### Save features ####

out_path <- file.path(output_dir, "Meta_features.csv")
write.csv(meta_features, out_path, row.names = FALSE)
cat("Wrote", nrow(meta_features), "patient rows to", out_path, "\n")
