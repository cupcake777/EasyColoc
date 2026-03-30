#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

source("src/utils_runtime.R")

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

tmp_dir <- tempfile(pattern = "easycoloc_output_tools_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "abf"), showWarnings = FALSE)
dir.create(file.path(tmp_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(tmp_dir, "rds"), showWarnings = FALSE)
dir.create(file.path(tmp_dir, "susie"), showWarnings = FALSE)

log_file <- file.path(tmp_dir, "run_easycoloc_full_20990101_000000.log")
writeLines(
  c(
    "Processing GWAS Dataset: SMOKE_GWAS",
    "[SUM] Merge complete!",
    "[REPORT] Interactive report saved: smoke_report.html"
  ),
  log_file
)

initialize_runtime_tracker(tmp_dir, enabled = TRUE, config_fingerprint = "smoke_fp")
register_active_run(log_file = log_file, run_label = "output_tools_smoke")
write_runtime_heartbeat(stage = "smoke", message_text = "output-dir tools smoke test")
runtime_update_task(
  task_id = build_task_id("SMOKE_GWAS", "chr1:1000:rs1", "prenatal", "GENE1"),
  status = "completed",
  gwas_id = "SMOKE_GWAS",
  locus_id = "chr1:1000:rs1",
  qtl_id = "prenatal",
  phenotype = "GENE1",
  result_file = file.path(tmp_dir, "abf", "SMOKE_GWAS_rs1_locus_results.csv"),
  pp4 = 0.87,
  n_snps = 31
)

fwrite(
  data.table(
    GWAS_ID = "SMOKE_GWAS",
    QTL_ID = "prenatal",
    Locus = "rs1",
    Phenotype = "GENE1",
    PP4 = 0.87,
    n_snps = 31
  ),
  file.path(tmp_dir, "abf", "SMOKE_GWAS_rs1_locus_results.csv")
)

Sys.setenv(
  EASYCOLOC_GLOBAL_CONFIG = file.path(tmp_dir, "missing_global.yaml"),
  EASYCOLOC_GWAS_CONFIG = file.path(tmp_dir, "missing_gwas.yaml"),
  EASYCOLOC_QTL_CONFIG = file.path(tmp_dir, "missing_qtl.yaml")
)
on.exit(
  Sys.unsetenv(c("EASYCOLOC_GLOBAL_CONFIG", "EASYCOLOC_GWAS_CONFIG", "EASYCOLOC_QTL_CONFIG")),
  add = TRUE
)

snapshot_file <- tempfile(pattern = "easycoloc_snapshot_", fileext = ".json")
monitor_out <- system2(
  "Rscript",
  c("tools/monitor_easycoloc.R", "--snapshot-file", snapshot_file, tmp_dir),
  stdout = TRUE,
  stderr = TRUE
)
assert_true(any(grepl("\\[MONITOR\\] output_dir:", monitor_out)), "monitor output missing output_dir")
assert_true(file.exists(snapshot_file), "monitor snapshot override file missing")

status_out <- system2(
  "Rscript",
  c("tools/summarize_run_status.R", tmp_dir),
  stdout = TRUE,
  stderr = TRUE
)
assert_true(any(grepl("SMOKE_GWAS", status_out)), "status output missing inferred GWAS id")
assert_true(file.exists(file.path(tmp_dir, "run_status_summary.txt")), "status text summary missing")

manifest_out <- system2(
  "Rscript",
  c("tools/build_output_manifest.R", tmp_dir),
  stdout = TRUE,
  stderr = TRUE
)
assert_true(any(grepl("\\[MANIFEST\\] total files:", manifest_out)), "manifest output missing summary")
assert_true(file.exists(file.path(tmp_dir, "output_manifest.tsv")), "manifest file missing")

check_out <- suppressWarnings(system2(
  "Rscript",
  c("tools/check_run_completion.R", tmp_dir),
  stdout = TRUE,
  stderr = TRUE
))
check_status <- attr(check_out, "status")
if (is.null(check_status)) {
  check_status <- 0L
}
assert_true(check_status == 1L, "completion check should report incomplete status for smoke output")
assert_true(any(grepl("\\[CHECK\\] status: INCOMPLETE", check_out)), "completion check output missing status")

cat("[SMOKE] output-dir tools smoke test passed\n")
cat("[SMOKE] output dir:", tmp_dir, "\n")
