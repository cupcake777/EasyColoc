#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
source("src/utils_config.R")

cfg_bundle <- easycoloc_read_configs()
output_dir <- normalizePath(cfg_bundle$global$output_dir, mustWork = FALSE)
requested_log <- NULL

if (length(args) >= 1 && nzchar(args[1])) {
  output_dir <- normalizePath(args[1], mustWork = FALSE)
}

if (length(args) >= 2 && nzchar(args[2])) {
  requested_log <- normalizePath(args[2], mustWork = FALSE)
}

latest_log <- NULL
if (!is.null(requested_log) && file.exists(requested_log)) {
  latest_log <- requested_log
} else if (dir.exists(output_dir)) {
  timestamped_logs <- list.files(
    output_dir,
    pattern = "^run_easycoloc_full_[0-9]{8}_[0-9]{4,6}\\.log$",
    full.names = TRUE
  )
  fallback_logs <- list.files(
    output_dir,
    pattern = "^run_easycoloc_full.*\\.log$",
    full.names = TRUE
  )

  log_candidates <- if (length(timestamped_logs) > 0) timestamped_logs else fallback_logs
  if (length(log_candidates) > 0) {
    log_info <- file.info(log_candidates)
    latest_log <- rownames(log_info)[which.max(log_info$mtime)]
  }
}

log_lines <- character(0)
if (!is.null(latest_log) && file.exists(latest_log)) {
  log_lines <- readLines(latest_log, warn = FALSE)
}

has_complete_marker <- any(grepl("EasyColoc Analysis Complete!", log_lines, fixed = TRUE))
has_summary_marker <- any(grepl("[SUM] Merge complete!", log_lines, fixed = TRUE))
has_report_marker <- any(grepl("[REPORT] Interactive report saved:", log_lines, fixed = TRUE))

abf_count <- if (dir.exists(file.path(output_dir, "abf"))) {
  length(list.files(file.path(output_dir, "abf"), pattern = "_locus_results\\.csv$", full.names = TRUE))
} else 0L

rds_count <- if (dir.exists(file.path(output_dir, "rds"))) {
  length(list.files(file.path(output_dir, "rds"), pattern = "\\.rds$", full.names = TRUE))
} else 0L

plot_count <- if (dir.exists(file.path(output_dir, "plots"))) {
  length(list.files(file.path(output_dir, "plots"), pattern = "\\.(pdf|png)$", full.names = TRUE))
} else 0L

susie_count <- if (dir.exists(file.path(output_dir, "susie"))) {
  length(list.files(file.path(output_dir, "susie"), pattern = "_susie\\.csv$", full.names = TRUE))
} else 0L

manifest_file <- file.path(output_dir, "output_manifest.tsv")
manifest_summary_file <- file.path(output_dir, "output_manifest_summary.json")

status <- if (has_complete_marker && has_summary_marker && has_report_marker) {
  "COMPLETE"
} else if (length(log_lines) > 0) {
  "INCOMPLETE"
} else {
  "NOT_STARTED"
}

cat("[CHECK] status:", status, "\n")
cat("[CHECK] output_dir:", output_dir, "\n")
cat("[CHECK] latest_log:", ifelse(is.null(latest_log), "NA", latest_log), "\n")
cat("[CHECK] abf:", abf_count, "rds:", rds_count, "plots:", plot_count, "susie:", susie_count, "\n")
cat("[CHECK] manifest:", ifelse(file.exists(manifest_file), manifest_file, "missing"), "\n")
cat("[CHECK] manifest_summary:", ifelse(file.exists(manifest_summary_file), manifest_summary_file, "missing"), "\n")

if (length(log_lines) > 0) {
  cat("[CHECK] last_log_line:", tail(log_lines, 1), "\n")
}

quit(status = if (identical(status, "COMPLETE")) 0 else 1)
