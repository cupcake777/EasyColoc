#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("src/utils_format.R")

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

sumstats_dt <- data.table(
  CHR = c(1, 1, 1),
  POS = c(100L, 120L, 150L),
  P = c(1e-8, 2e-7, 3e-6)
)

lifted_positions <- data.table(
  V1 = c("chr1", "chr1"),
  V2 = c(1000L, 1300L),
  V3 = c(1001L, 1301L)
)

updated_dt <- apply_liftover_positions(
  sumstats_dt = sumstats_dt,
  chrom = "1",
  start_pos = 100L,
  end_pos = 150L,
  hg38_positions = lifted_positions
)

assert_true(identical(updated_dt$CHR, c("1", "1", "1")), "CHR should be updated to lifted chromosome")
assert_true(identical(updated_dt$POS, c(1001L, 1121L, 1301L)), "POS should be remapped into the lifted interval")
assert_true(identical(sumstats_dt$POS, c(100L, 120L, 150L)), "input data should remain unchanged")

cat("[SMOKE] harmonization fallback smoke test passed\n")
