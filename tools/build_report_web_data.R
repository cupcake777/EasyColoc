#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

read_flag <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx >= length(args)) {
    return(NULL)
  }
  args[[idx + 1L]]
}

results_dir <- read_flag("--results-dir")
project_name <- read_flag("--project-name")

if (is.null(results_dir) || !nzchar(results_dir)) {
  stop(
    "Usage: Rscript tools/build_report_web_data.R --results-dir PATH [--project-name NAME]",
    call. = FALSE
  )
}

source("src/utils_report_web.R")

out_file <- write_report_web_payload(
  results_dir = results_dir,
  project_name = project_name
)

message("[REPORT-WEB] payload: ", out_file)
