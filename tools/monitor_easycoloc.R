#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
source("src/utils_config.R")

cfg_bundle <- easycoloc_read_configs()
output_dir <- normalizePath(cfg_bundle$global$output_dir, mustWork = FALSE)
if (length(args) >= 1 && nzchar(args[1])) {
  output_dir <- normalizePath(args[1], mustWork = FALSE)
}

runtime_dir <- file.path(output_dir, "runtime")
active_run_file <- file.path(runtime_dir, "active_run.json")
heartbeat_file <- file.path(runtime_dir, "heartbeat.json")
task_state_file <- file.path(runtime_dir, "task_state.tsv")
snapshot_file <- file.path(runtime_dir, "monitor_snapshot.json")

read_json_safely <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

active_run <- read_json_safely(active_run_file)
heartbeat <- read_json_safely(heartbeat_file)
task_state <- if (file.exists(task_state_file)) {
  tryCatch(as.data.table(fread(task_state_file)), error = function(e) data.table())
} else {
  data.table()
}

log_candidates <- character(0)
if (dir.exists(output_dir)) {
  log_candidates <- list.files(
    output_dir,
    pattern = "^run_easycoloc_full.*\\.log$",
    full.names = TRUE
  )
}

active_log <- if (!is.null(active_run) && !is.null(active_run$log_file) && file.exists(active_run$log_file)) {
  active_run$log_file
} else if (length(log_candidates) > 0) {
  info <- file.info(log_candidates)
  rownames(info)[which.max(info$mtime)]
} else {
  NA_character_
}

last_log_line <- NA_character_
if (!is.na(active_log) && file.exists(active_log)) {
  log_lines <- readLines(active_log, warn = FALSE)
  if (length(log_lines) > 0) last_log_line <- tail(log_lines, 1)
}

count_outputs <- function(subdir, pattern) {
  dir_path <- file.path(output_dir, subdir)
  if (!dir.exists(dir_path)) return(0L)
  length(list.files(dir_path, pattern = pattern, full.names = TRUE))
}

task_counts <- if (nrow(task_state) > 0 && "status" %in% names(task_state)) {
  as.list(table(task_state$status, useNA = "ifany"))
} else {
  list()
}

snapshot <- list(
  output_dir = output_dir,
  active_log = active_log,
  active_run = active_run,
  heartbeat = heartbeat,
  last_log_line = last_log_line,
  outputs = list(
    abf = count_outputs("abf", "_locus_results\\.csv$"),
    rds = count_outputs("rds", "\\.rds$"),
    plots = count_outputs("plots", "\\.(pdf|png)$"),
    susie = count_outputs("susie", "_susie\\.csv$")
  ),
  task_counts = task_counts,
  task_rows = nrow(task_state)
)

if (dir.exists(runtime_dir)) {
  write_json(snapshot, snapshot_file, auto_unbox = TRUE, pretty = TRUE, null = "null")
}

cat("[MONITOR] output_dir:", output_dir, "\n")
cat("[MONITOR] active_log:", ifelse(is.na(active_log), "NA", active_log), "\n")
if (!is.null(active_run)) {
  cat("[MONITOR] active_run_pid:", ifelse(is.null(active_run$pid), "NA", active_run$pid), "\n")
  cat("[MONITOR] active_run_label:", ifelse(is.null(active_run$run_label), "NA", active_run$run_label), "\n")
}
if (!is.null(heartbeat)) {
  cat("[MONITOR] heartbeat_stage:", ifelse(is.null(heartbeat$stage), "NA", heartbeat$stage), "\n")
  cat("[MONITOR] heartbeat_time:", ifelse(is.null(heartbeat$timestamp), "NA", heartbeat$timestamp), "\n")
}
cat("[MONITOR] outputs:",
    sprintf("abf=%s rds=%s plots=%s susie=%s",
            snapshot$outputs$abf,
            snapshot$outputs$rds,
            snapshot$outputs$plots,
            snapshot$outputs$susie),
    "\n")
if (length(task_counts) > 0) {
  cat("[MONITOR] task_status:",
      paste(sprintf("%s=%s", names(task_counts), unlist(task_counts)), collapse = " "),
      "\n")
}
cat("[MONITOR] last_log_line:", ifelse(is.na(last_log_line), "NA", last_log_line), "\n")
