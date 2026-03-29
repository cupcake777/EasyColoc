#!/usr/bin/env Rscript

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

tmp_dir <- tempfile(pattern = "easycoloc_demo_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
project_dir <- file.path(tmp_dir, "demo_project")

cmd <- c("easycoloc", "bootstrap-refs", "--demo", project_dir, "--force", "--run")
out <- system2("bash", cmd, stdout = TRUE, stderr = TRUE)

assert_true(file.exists(file.path(project_dir, "config", "global.yaml")), "demo global config missing")
assert_true(file.exists(file.path(project_dir, "data", "gwas", "demo_gwas.tsv.gz")), "demo GWAS file missing")
assert_true(file.exists(file.path(project_dir, "data", "qtl", "all_pairs_demo.txt.gz")), "demo allPairs file missing")
assert_true(file.exists(file.path(project_dir, "data", "qtl", "all_pairs_demo.txt.gz.tbi")), "demo allPairs tabix missing")
assert_true(file.exists(file.path(project_dir, "refs", "plink", "hg38", "toy_chr22.bed")), "demo PLINK bed missing")
assert_true(file.exists(file.path(project_dir, "refs", "plink", "keep", "DEMO.sample")), "demo keep file missing")
assert_true(file.exists(file.path(project_dir, "results", "coloc_report.html")), "demo HTML report missing")
assert_true(any(grepl("demo_project", out, fixed = TRUE)), "demo bootstrap output missing project path")
cat("[SMOKE] bootstrap demo smoke test passed\n")
