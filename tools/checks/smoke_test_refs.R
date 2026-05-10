#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

tmp_dir <- tempfile(pattern = "easycoloc_refs_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "config"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "data", "gwas"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "data", "qtl"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "ref"), recursive = TRUE, showWarnings = FALSE)

plink_prefix <- file.path(tmp_dir, "ref", "toy_plink")
invisible(file.create(paste0(plink_prefix, c(".bed", ".bim", ".fam"))))
invisible(file.create(
  file.path(
    tmp_dir,
    "ref",
    c("keep.sample", "gene.gtf.gz", "ref19.fa.gz", "ref38.fa.gz", "dbsnp19.gz", "dbsnp38.gz", "hg19ToHg38.over.chain.gz")
  )
))
dir.create(file.path(tmp_dir, "ref", "af_ref"), recursive = TRUE, showWarnings = FALSE)
invisible(file.create(file.path(tmp_dir, "ref", "af_ref", "1KG_hg19_EUR_AF.tsv.gz")))
invisible(file.create(file.path(tmp_dir, "data", "gwas", "example.tsv.gz")))
invisible(file.create(file.path(tmp_dir, "data", "qtl", c("all_pairs.txt.gz", "sig_pairs.txt.gz"))))
writeLines(
  c(
    "Type,NumberRNASeqandGTSamples,allPairsTabixFilename,sigPairsTabixFilename",
    "example_eqtl,500,all_pairs.txt.gz,sig_pairs.txt.gz"
  ),
  file.path(tmp_dir, "data", "qtl", "QTL_summary.csv")
)

global_cfg <- list(
  project_root = "..",
  output_dir = "results",
  temp_dir = "temp",
  plink_hg38 = "ref/toy_plink",
  plink_keep = "ref/keep.sample",
  hash_table_dir = "",
  recom = "",
  gene_anno = "ref/gene.gtf.gz",
  ref_genome_hg19 = "ref/ref19.fa.gz",
  ref_genome_hg38 = "ref/ref38.fa.gz",
  `1kg_af` = "ref/af_ref",
  dbsnp_hg19 = "ref/dbsnp19.gz",
  dbsnp_hg38 = "ref/dbsnp38.gz",
  harmonize_dir = "harmony",
  clump = list(p1 = 5e-8, p2 = 5e-8, r2 = 0.1, kb = 1000),
  coloc_settings = list(pp4_threshold = 0.7, susie_threshold = 0.5, min_snps = 10),
  harmonization_settings = list(liftover_chain = "ref/hg19ToHg38.over.chain.gz")
)
yaml::write_yaml(global_cfg, file.path(tmp_dir, "config", "global.yml"))

gwas_cfg <- list(
  datasets = list(
    list(
      id = "EXAMPLE",
      name = "Example",
      file = "data/gwas/example.tsv.gz",
      build = "hg19",
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
yaml::write_yaml(gwas_cfg, file.path(tmp_dir, "config", "gwas.yml"))

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
yaml::write_yaml(qtl_cfg, file.path(tmp_dir, "config", "qtl.yml"))

cmd <- c(
  "tools/list_reference_requirements.R",
  "--global", file.path(tmp_dir, "config", "global.yml"),
  "--gwas", file.path(tmp_dir, "config", "gwas.yml"),
  "--qtl", file.path(tmp_dir, "config", "qtl.yml"),
  "--include-qtl-files"
)
refs_out <- system2("Rscript", cmd, stdout = TRUE, stderr = TRUE)

assert_true(any(grepl("\\[REF\\].*plink_hg38", refs_out)), "plink_hg38 requirement missing")
assert_true(any(grepl("\\[REF\\].*qtl_summary", refs_out)), "qtl_summary requirement missing")
assert_true(any(grepl("qtl_allpairs_example_eqtl", refs_out, fixed = TRUE)), "resolved QTL allPairs entry missing")
cat("[SMOKE] reference listing smoke test passed\n")
