#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
source("src/utils_config.R")
cfg_bundle <- easycoloc_read_configs()
cfg_global <- cfg_bundle$global
cfg_gwas <- cfg_bundle$gwas

output_dir <- normalizePath(cfg_global$output_dir, mustWork = FALSE)
if (length(args) >= 1 && nzchar(args[1])) {
  output_dir <- normalizePath(args[1], mustWork = FALSE)
}

timestamped_logs <- if (dir.exists(output_dir)) {
  list.files(output_dir, pattern = "^run_easycoloc_full_[0-9]{8}_[0-9]{4,6}\\.log$", full.names = TRUE)
} else {
  character(0)
}

latest_log <- NULL
if (length(timestamped_logs) > 0) {
  info <- file.info(timestamped_logs)
  latest_log <- rownames(info)[which.max(info$mtime)]
}

runtime_dir <- file.path(output_dir, "runtime")
active_run_file <- file.path(runtime_dir, "active_run.json")
heartbeat_file <- file.path(runtime_dir, "heartbeat.json")
task_state_file <- file.path(runtime_dir, "task_state.tsv")

read_json_safely <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

active_run <- read_json_safely(active_run_file)
heartbeat <- read_json_safely(heartbeat_file)
if (!is.null(active_run) && !is.null(active_run$log_file) && file.exists(active_run$log_file)) {
  latest_log <- active_run$log_file
}

log_lines <- if (!is.null(latest_log) && file.exists(latest_log)) readLines(latest_log, warn = FALSE) else character(0)
last_started_gwas <- {
  hits <- grep("^Processing GWAS Dataset:", log_lines, value = TRUE)
  if (length(hits) > 0) sub("^Processing GWAS Dataset: ", "", tail(hits, 1)) else NA_character_
}

gwas_ids <- vapply(cfg_gwas$datasets, `[[`, character(1), "id")

count_prefix_files <- function(subdir, pattern_suffix) {
  dir_path <- file.path(output_dir, subdir)
  if (!dir.exists(dir_path)) return(setNames(integer(length(gwas_ids)), gwas_ids))
  files <- list.files(dir_path, full.names = FALSE)
  setNames(
    vapply(gwas_ids, function(id) {
      sum(grepl(paste0("^", id, ".*", pattern_suffix, "$"), files))
    }, integer(1)),
    gwas_ids
  )
}

abf_counts <- count_prefix_files("abf", "_locus_results\\.csv")
rds_counts <- count_prefix_files("rds", "\\.rds")
plot_counts <- count_prefix_files("plots", "\\.(pdf|png)")
susie_counts <- count_prefix_files("susie", "_susie\\.csv")

summary_dt <- data.table(
  gwas_id = gwas_ids,
  abf = as.integer(abf_counts[gwas_ids]),
  rds = as.integer(rds_counts[gwas_ids]),
  plots = as.integer(plot_counts[gwas_ids]),
  susie = as.integer(susie_counts[gwas_ids])
)
summary_dt[, started := gwas_id %in% gsub("^.*: ", "", grep("^Processing GWAS Dataset:", log_lines, value = TRUE))]
summary_dt[, current := !is.na(last_started_gwas) & gwas_id == last_started_gwas]

task_status_dt <- if (file.exists(task_state_file)) {
  tryCatch(as.data.table(fread(task_state_file)), error = function(e) data.table())
} else {
  data.table()
}

txt_file <- file.path(output_dir, "run_status_summary.txt")
json_file <- file.path(output_dir, "run_status_summary.json")

writeLines(
  c(
    paste0("output_dir: ", output_dir),
    paste0("latest_log: ", ifelse(is.null(latest_log), "NA", latest_log)),
    paste0("last_started_gwas: ", ifelse(is.na(last_started_gwas), "NA", last_started_gwas)),
    paste0("active_run_pid: ", ifelse(is.null(active_run$pid), "NA", active_run$pid)),
    paste0("heartbeat_stage: ", ifelse(is.null(heartbeat$stage), "NA", heartbeat$stage)),
    paste0("heartbeat_time: ", ifelse(is.null(heartbeat$timestamp), "NA", heartbeat$timestamp)),
    "",
    capture.output(print(summary_dt)),
    "",
    if (nrow(task_status_dt) > 0) {
      c("task_status_counts:", capture.output(print(task_status_dt[, .N, by = status][order(status)])))
    } else {
      "task_status_counts: none"
    }
  ),
  txt_file
)

write_json(
  list(
    output_dir = output_dir,
    latest_log = latest_log,
    last_started_gwas = last_started_gwas,
    active_run = active_run,
    heartbeat = heartbeat,
    summary = summary_dt,
    task_status_counts = if (nrow(task_status_dt) > 0) task_status_dt[, .N, by = status][order(status)] else data.table()
  ),
  json_file,
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("[STATUS] wrote:", txt_file, "\n")
cat("[STATUS] wrote:", json_file, "\n")
print(summary_dt)
