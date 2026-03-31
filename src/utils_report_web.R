suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

resolve_current_script_dir <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else as.character(ofile)
    },
    character(1)
  )
  frame_files <- frame_files[nzchar(frame_files)]

  if (length(frame_files) == 0) {
    return(getwd())
  }

  dirname(normalizePath(frame_files[[length(frame_files)]], mustWork = FALSE))
}

utils_report_web_dir <- resolve_current_script_dir()
source(file.path(utils_report_web_dir, "utils_report.R"))

read_json_if_present <- function(path) {
  if (!file.exists(path)) {
    return(list(missing = TRUE, path = basename(path)))
  }

  parsed <- tryCatch(
    read_json(path, simplifyVector = FALSE),
    error = function(e) list(invalid_json = TRUE, path = basename(path), error = e$message)
  )
  parsed
}

normalize_report_results <- function(results_dt) {
  if (is.null(results_dt) || nrow(results_dt) == 0) {
    return(list())
  }

  rows <- split(results_dt, seq_len(nrow(results_dt)))
  out <- vector("list", length(rows))

  for (i in seq_along(rows)) {
    row <- rows[[i]]
    out[[i]] <- list(
      gwas_id = as.character(report_value(row, "GWAS_ID", default = "")),
      qtl_id = as.character(report_value(row, "QTL_ID", default = "")),
      locus = as.character(report_value(row, "Locus", default = "")),
      phenotype = as.character(report_value(row, "Gene", fallback = "Phenotype", default = "")),
      pp4 = as.numeric(report_value(row, "PP4", default = NA_real_)),
      n_snps = as.integer(report_value(row, "n_snps", default = NA_integer_)),
      source_file = as.character(report_value(row, "source_file", default = ""))
    )
  }

  out
}

index_report_assets <- function(results_dir) {
  plot_dir <- file.path(results_dir, "plots")
  plot_files <- if (dir.exists(plot_dir)) list.files(plot_dir, full.names = FALSE) else character(0)
  plots <- lapply(plot_files, function(name) {
    list(
      name = name,
      rel_path = file.path("plots", name)
    )
  })

  top_level_entries <- list.files(results_dir, full.names = FALSE, recursive = FALSE)
  top_level_files <- top_level_entries[!dir.exists(file.path(results_dir, top_level_entries))]
  download_candidates <- top_level_files[grepl("\\.(csv|tsv|txt|json)$", top_level_files, ignore.case = TRUE)]
  download_files <- sort(setdiff(download_candidates, c("heartbeat.json", "monitor_snapshot.json")))

  downloads <- lapply(download_files, function(name) {
    list(
      name = name,
      rel_path = name
    )
  })

  list(
    plots = plots,
    downloads = downloads
  )
}

build_report_web_payload <- function(results_dir, project_name = NULL) {
  abf_files <- list.files(
    file.path(results_dir, "abf"),
    pattern = "_results\\.csv$",
    full.names = TRUE
  )

  results_dt <- merge_results(abf_files)
  if (is.null(results_dt) || nrow(results_dt) == 0) {
    stop("No ABF result files found under ", file.path(results_dir, "abf"), call. = FALSE)
  }

  project_name_resolved <- if (is.null(project_name) || !nzchar(project_name)) {
    basename(normalizePath(results_dir, mustWork = FALSE))
  } else {
    project_name
  }

  list(
    meta = list(
      project_name = project_name_resolved,
      results_dir = normalizePath(results_dir, mustWork = FALSE),
      generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      report_version = "report-web-v1"
    ),
    summary = get_summary_stats(results_dt),
    results = normalize_report_results(results_dt),
    assets = index_report_assets(results_dir),
    runtime = list(
      heartbeat = read_json_if_present(file.path(results_dir, "heartbeat.json")),
      monitor_snapshot = read_json_if_present(file.path(results_dir, "monitor_snapshot.json"))
    ),
    warnings = list()
  )
}

write_report_web_payload <- function(results_dir, project_name = NULL) {
  output_dir <- file.path(results_dir, "report_web")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  payload <- build_report_web_payload(results_dir, project_name = project_name)
  out_file <- file.path(output_dir, "report-data.json")
  write_json(payload, out_file, auto_unbox = TRUE, pretty = TRUE, null = "null")

  out_file
}
