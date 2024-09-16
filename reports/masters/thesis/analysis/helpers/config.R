library(R6)

Config <- R6Class("Config",  list(
  experiment_path = NULL,
  samplesheet_path = NULL,
  sample_depth = NA,
  repetition = NA,

  initialize = function(experiment_path, samplesheet_path, sample_depth, repetition) {
    stopifnot(is.character(experiment_path), length(experiment_path) == 1)
    stopifnot(is.character(samplesheet_path), length(samplesheet_path) == 1)
    stopifnot(is.numeric(sample_depth), length(sample_depth) == 1)
    stopifnot(is.numeric(repetition), length(repetition) == 1)

    self$experiment_path <- experiment_path
    self$samplesheet_path <- samplesheet_path
    self$sample_depth <- sample_depth
    self$repetition <- repetition
  }
))

config <- Config$new(
  # experiment_path = "../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14",
  # experiment_path = "../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14-NC",
  experiment_path = "../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-09-12",
  samplesheet_path = "../../../../experiments/66-fungal-isolate-ONT/samplesheet-unite-tax.csv",
  sample_depth = 2000,
  repetition = 2
)

read_samplesheet <- function(config) {
    samplesheet <- read.csv(config$samplesheet_path, header = TRUE, row.names = 1)
    rownames(samplesheet) <- paste0("barcode", rownames(samplesheet))
    samplesheet
}
