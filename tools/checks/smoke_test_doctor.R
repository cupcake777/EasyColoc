#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

source("src/utils_config.R")
source("tools/checks/smoke_project_fixture.R")

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

env_root <- tempfile(pattern = "easycoloc_env_root_")
Sys.setenv(EASYCOLOC_SMOKE_ROOT = env_root)
expanded_paths <- easycoloc_expand_env_vars(c(
  "$EASYCOLOC_SMOKE_ROOT/ref",
  "${EASYCOLOC_SMOKE_ROOT}/data",
  "plain/path"
))
assert_true(identical(expanded_paths[[1]], file.path(env_root, "ref")), "bare env var expansion failed")
assert_true(identical(expanded_paths[[2]], file.path(env_root, "data")), "braced env var expansion failed")
assert_true(identical(expanded_paths[[3]], "plain/path"), "plain path should not be rewritten")

fixture <- create_smoke_project_fixture(
  pattern = "easycoloc_doctor_",
  gwas_build = "hg38",
  recom = "ref/recom_map",
  include_results = TRUE,
  qtl_paths = "absolute"
)

cmd <- c(
  "tools/doctor_easycoloc.R",
  "--global", fixture$global,
  "--gwas", fixture$gwas,
  "--qtl", fixture$qtl
)
doctor_out <- system2("Rscript", cmd, stdout = TRUE, stderr = TRUE)

assert_true(any(grepl("\\[DOCTOR\\] summary:", doctor_out)), "doctor summary line missing")
assert_true(!any(grepl("\\[FAIL\\]", doctor_out)), "doctor unexpectedly reported FAIL")
cat("[SMOKE] doctor smoke test passed\n")
