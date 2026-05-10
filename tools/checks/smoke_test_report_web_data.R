#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

script_arg <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
repo_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE)
build_script <- file.path(repo_root, "tools", "build_report_web_data.R")

tmp_dir <- tempfile(pattern = "easycoloc_report_web_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "abf"), showWarnings = FALSE)
dir.create(file.path(tmp_dir, "susie"), showWarnings = FALSE)
dir.create(file.path(tmp_dir, "plots"), showWarnings = FALSE)

fwrite(
  data.table(
    GWAS_ID = "GWAS_A",
    QTL_ID = "postnatal",
    Locus = "rs100",
    Phenotype = "GENE1",
    PP4 = 0.92,
    n_snps = 41
  ),
  file.path(tmp_dir, "abf", "GWAS_A_rs100_locus_results.csv")
)

fwrite(
  data.table(
    cs_id = 1L,
    variable = "rs100",
    variable_prob = 0.78,
    cs_log10bf = 11.2
  ),
  file.path(tmp_dir, "susie", "GWAS_A_postnatal_GENE1_susie.csv")
)

writeLines("plot placeholder", file.path(tmp_dir, "plots", "GWAS_A_rs100_panel.png"))
writeLines(
  '{"stage":["pipeline_complete"],"message_text":"done","codes":[101]}',
  file.path(tmp_dir, "heartbeat.json")
)
fwrite(data.table(metric = "ok", value = 1), file.path(tmp_dir, "extra_download.tsv"), sep = "\t")

run_dir <- tempfile(pattern = "easycoloc_report_web_run_")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(run_dir)

status <- system2(
  "Rscript",
  c(build_script, "--results-dir", tmp_dir),
  stdout = TRUE,
  stderr = TRUE
)

payload_path <- file.path(tmp_dir, "report_web", "report-data.json")
assert_true(file.exists(payload_path), "report-data.json was not created")

payload <- read_json(payload_path, simplifyVector = FALSE)
assert_true(identical(payload$summary$total_tests, 1L), "summary$total_tests mismatch")
assert_true(is.null(names(payload$results)), "results must be an array, not a named object")
assert_true(length(payload$results) == 1L, "results length mismatch")
assert_true(is.null(names(payload$assets$plots)), "assets$plots must be an array")
assert_true(length(payload$assets$plots) == 1L, "plot index mismatch")
assert_true(is.null(names(payload$assets$downloads)), "assets$downloads must be an array")
download_names <- vapply(payload$assets$downloads, `[[`, "", "name")
assert_true("extra_download.tsv" %in% download_names, "extra download file missing from assets$downloads")
assert_true(is.list(payload$runtime$heartbeat$stage), "heartbeat stage should remain an array")
assert_true(length(payload$runtime$heartbeat$stage) == 1L, "heartbeat stage array length mismatch")
assert_true(identical(payload$runtime$heartbeat$stage[[1]], "pipeline_complete"), "heartbeat stage mismatch")
assert_true(is.list(payload$runtime$heartbeat$codes), "heartbeat codes should remain an array")
assert_true(length(payload$runtime$heartbeat$codes) == 1L, "heartbeat codes array length mismatch")
assert_true(isTRUE(all.equal(as.numeric(payload$runtime$heartbeat$codes[[1]]), 101)), "heartbeat codes value mismatch")

cat("[SMOKE] report-web data smoke test passed\n")
