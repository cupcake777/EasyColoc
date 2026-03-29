#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
})

utils_files <- list.files("src", pattern = "^utils_.*\\.R$", full.names = TRUE)
invisible(lapply(utils_files, source))

cfg_global <- yaml::read_yaml("config/global.yaml")

make_synthetic_locus <- function() {
  set.seed(20260327)

  n_bg <- 80
  bg_pos <- sort(sample(seq(43120000, 43580000, by = 500), n_bg))
  bg_p <- runif(n_bg, min = 5e-3, max = 0.08)

  signal_pos <- c(43562000, 43564000, 43566000, 43568000, 43570000, 43572000)
  signal_p <- c(2e-6, 8e-7, 2e-8, 7e-8, 3e-6, 9e-6)
  signal_rs <- c("rs2819340", "rs3791138", "rs10890255", "rs3791137", "rs2004899", "rs2819341")

  df <- data.frame(
    rsid = c(
      paste0("rsBG", seq_len(n_bg)),
      signal_rs
    ),
    CHR.qtl = "chr1",
    POS.qtl = c(bg_pos, signal_pos),
    P.qtl = c(bg_p, signal_p),
    P.gwas = c(runif(n_bg, min = 1e-3, max = 0.2), c(3e-5, 6e-6, 1e-9, 5e-6, 9e-5, 7e-5)),
    stringsAsFactors = FALSE
  )

  credible_set <- data.frame(
    snp = c("rs2819340", "rs3791138", "rs2004899"),
    SNP.PP.H4 = c(0.42, 0.31, 0.18),
    stringsAsFactors = FALSE
  )

  list(
    merged_data = df,
    credible_set = credible_set,
    lead_snp = "rs10890255",
    phenotype_info = "ENST00000654683.1|CCDC30|chr1:42482663-42484158|+",
    qtl_type = "postnatal"
  )
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

example_data <- make_synthetic_locus()

feature_region <- parse_feature_region(example_data$phenotype_info)
assert_true(!is.null(feature_region), "parse_feature_region() returned NULL")
assert_true(feature_region$chrom == "chr1", "feature_region chrom parse failed")
assert_true(feature_region$start == 42482663, "feature_region start parse failed")
assert_true(feature_region$end == 42484158, "feature_region end parse failed")

window_info <- resolve_plot_window(
  leadSNP_DF = example_data$merged_data,
  pos_col = "POS.qtl",
  chr_col = "CHR.qtl",
  lead_snp_val = example_data$lead_snp,
  plot_window_bp = 200000,
  phenotype_info = example_data$phenotype_info
)
assert_true(!is.null(window_info), "resolve_plot_window() returned NULL")
assert_true(window_info$chrom_label == "chr1", "resolve_plot_window() chrom label failed")
assert_true(isTRUE(window_info$expanded), "resolve_plot_window() should expand to cover phenotype transcript")
assert_true(window_info$min_pos == 42482663, "resolve_plot_window() expanded min_pos is incorrect")
assert_true(window_info$max_pos == 43766000, "resolve_plot_window() expanded max_pos is incorrect")

assign("lead_SNP", example_data$lead_snp, envir = .GlobalEnv)
assign("geneSymbol", example_data$phenotype_info, envir = .GlobalEnv)
assign("plink_bfile", NULL, envir = .GlobalEnv)
assign("trait", "SMOKE_GWAS", envir = .GlobalEnv)

plot_obj <- plot_qtl_association(
  qtl_all_chrom = "CHR.qtl",
  qtl_all_pvalue = "P.qtl",
  leadSNP_DF = example_data$merged_data,
  ld_df = NULL,
  gtf_path = cfg_global$gene_anno,
  region_recomb = NULL,
  recomb_path = cfg_global$recom,
  show_lead_line = FALSE,
  qtl_type = example_data$qtl_type,
  phenotype_info = example_data$phenotype_info,
  significance_threshold = 5e-8,
  significance_label = NULL,
  title_phenotype_field = "gene",
  plot_width = 6.6,
  plot_height = 4.9,
  credible_set = example_data$credible_set,
  plot_window_bp = 200000
)

tmp_pdf <- tempfile(fileext = ".pdf")
tmp_png <- tempfile(fileext = ".png")
save_res <- save_plot_with_fallback(
  plot_obj = plot_obj,
  pdf_path = tmp_pdf,
  png_path = tmp_png,
  plot_width = 6.6,
  plot_height = 4.9,
  png_scale = 1,
  png_dpi = 200
)

output_path <- save_res$path
assert_true(file.exists(output_path), "smoke plot output file was not created")
assert_true(file.info(output_path)$size > 5000, "smoke plot output file is suspiciously small")

cat("[SMOKE] plotting smoke test passed\n")
cat("[SMOKE] output:", output_path, "\n")
