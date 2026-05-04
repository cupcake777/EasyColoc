#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

source("src/utils_config.R")

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

tmp_dir <- tempfile(pattern = "easycoloc_doctor_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "config"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "data", "gwas"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "data", "qtl"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "ref"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "results"), recursive = TRUE, showWarnings = FALSE)

plink_prefix <- file.path(tmp_dir, "ref", "toy_plink")
invisible(file.create(paste0(plink_prefix, c(".bed", ".bim", ".fam"))))
invisible(file.create(
  file.path(
    tmp_dir,
    "ref",
    c(
      "keep.sample", "recom_map", "gene.gtf.gz", "ref19.fa.gz", "ref38.fa.gz",
      "af_ref", "dbsnp19.gz", "dbsnp38.gz", "hg19ToHg38.over.chain.gz"
    )
  )
))
invisible(file.create(file.path(tmp_dir, "data", "gwas", "example.tsv.gz")))
writeLines(
  c(
    "Type,NumberRNASeqandGTSamples,allPairsTabixFilename,sigPairsTabixFilename",
    paste("example_eqtl,500", file.path(tmp_dir, "data", "qtl", "all_pairs.txt.gz"), file.path(tmp_dir, "data", "qtl", "sig_pairs.txt.gz"), sep = ",")
  ),
  file.path(tmp_dir, "data", "qtl", "QTL_summary.csv")
)
invisible(file.create(file.path(tmp_dir, "data", "qtl", c("all_pairs.txt.gz", "sig_pairs.txt.gz"))))

global_cfg <- list(
  project_root = "..",
  output_dir = "results",
  temp_dir = "temp",
  plink_hg38 = "ref/toy_plink",
  plink_keep = "ref/keep.sample",
  hash_table_dir = "",
  recom = "ref/recom_map",
  gene_anno = "ref/gene.gtf.gz",
  ref_genome_hg19 = "ref/ref19.fa.gz",
  ref_genome_hg38 = "ref/ref38.fa.gz",
  `1kg_af` = "ref/af_ref",
  dbsnp_hg19 = "ref/dbsnp19.gz",
  dbsnp_hg38 = "ref/dbsnp38.gz",
  harmonize_dir = "harmony",
  clump = list(p1 = 5e-8, p2 = 5e-8, r2 = 0.1, kb = 1000),
  coloc_settings = list(pp4_threshold = 0.7, susie_threshold = 0.5, min_snps = 10),
  harmonization_settings = list(env_name = "gwaslab", liftover_chain = "ref/hg19ToHg38.over.chain.gz")
)
yaml::write_yaml(global_cfg, file.path(tmp_dir, "config", "global.yaml"))

gwas_cfg <- list(
  datasets = list(
    list(
      id = "EXAMPLE",
      name = "Example",
      file = "data/gwas/example.tsv.gz",
      build = "hg38",
      pop = "EUR",
      type = "cc",
      sample_size_n = 10000,
      prop = 0.2,
      columns = list(
        snp = "SNP",
        chrom = "CHR",
        pos = "BP",
        a1 = "A1",
        a2 = "A2",
        pval = "P",
        beta = "BETA",
        se = "SE"
      )
    )
  )
)
yaml::write_yaml(gwas_cfg, file.path(tmp_dir, "config", "gwas.yaml"))

qtl_cfg <- list(
  qtl_info = list(
    build = "hg38",
    file = "data/qtl/QTL_summary.csv",
    columns = list(
      id = "Type",
      sample_size = "NumberRNASeqandGTSamples",
      all_filename = "allPairsTabixFilename",
      sig_filename = "sigPairsTabixFilename"
    )
  ),
  QTL_all_header = c("chr", "pos", "phenotype_id", "variant_id", "pval_nominal", "slope", "slope_se"),
  QTL_cols = list(
    phenotype = "phenotype_id",
    chrom = "chr",
    pos = "pos",
    pval = "pval_nominal",
    beta = "slope",
    se = "slope_se"
  )
)
yaml::write_yaml(qtl_cfg, file.path(tmp_dir, "config", "qtl.yaml"))

cmd <- c(
  "tools/doctor_easycoloc.R",
  "--global", file.path(tmp_dir, "config", "global.yaml"),
  "--gwas", file.path(tmp_dir, "config", "gwas.yaml"),
  "--qtl", file.path(tmp_dir, "config", "qtl.yaml")
)
doctor_out <- system2("Rscript", cmd, stdout = TRUE, stderr = TRUE)

assert_true(any(grepl("\\[DOCTOR\\] summary:", doctor_out)), "doctor summary line missing")
assert_true(!any(grepl("\\[FAIL\\]", doctor_out)), "doctor unexpectedly reported FAIL")
cat("[SMOKE] doctor smoke test passed\n")
