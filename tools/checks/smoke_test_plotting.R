#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
})

source("src/utils_helpers.R")
source("src/utils_plot.R")
source("src/utils_output.R")

cfg_global <- yaml::read_yaml("config/global.yml")
fixture_dir <- file.path("examples", "fixtures")

resolve_existing_file <- function(paths) {
  paths <- Filter(function(path) !is.null(path) && nzchar(path), paths)
  for (path in paths) {
    if (file.exists(path)) return(path)
  }
  NULL
}

resolve_recomb_prefix <- function(prefixes, chrom = "1") {
  prefixes <- Filter(function(path) !is.null(path) && nzchar(path), prefixes)
  for (prefix in prefixes) {
    txt_path <- paste0(prefix, "_recombination_map_hapmap_format_hg38_chr_", chrom, ".txt")
    bed_path <- paste0(prefix, "_recombination_map_hg38_chr_", chrom, ".bed")
    if (file.exists(txt_path) || file.exists(bed_path)) return(prefix)
  }
  NULL
}

gtf_path <- resolve_existing_file(c(
  file.path(fixture_dir, "annotation", "smoke_hg38_chr1.gtf"),
  cfg_global$gene_anno
))
recomb_path <- resolve_recomb_prefix(c(
  file.path(fixture_dir, "recomb", "hg38", "CHB", "CHB"),
  cfg_global$recom
))

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

make_dense_corner_points <- function() {
  data.frame(
    x = c(seq(0.02, 0.24, length.out = 28), seq(0.68, 0.86, length.out = 4)),
    y = c(seq(0.66, 0.98, length.out = 28), seq(0.12, 0.22, length.out = 4)),
    is_lead = c(TRUE, rep(FALSE, 31)),
    is_credible = c(rep(TRUE, 4), rep(FALSE, 28)),
    stringsAsFactors = FALSE
  )
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

assert_true(!is.null(gtf_path), "smoke plotting test requires a resolvable gene annotation fixture")
assert_true(!is.null(recomb_path), "smoke plotting test requires a resolvable CHB recombination fixture")

example_data <- make_synthetic_locus()
dense_corner_points <- make_dense_corner_points()

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

legend_layout <- choose_internal_legend_position(
  points_df = dense_corner_points,
  xlim = c(0, 1),
  ylim = c(0, 1),
  legend_nrow = 6,
  legend_labels = c("Lead SNP", "95% CS", "0.8-1.0", "0.6-0.8", "0.4-0.6", "< 0.2"),
  plot_width = 6.6,
  plot_height = 4.9
)
assert_true(
  !(legend_layout$position[1] < 0.5 && legend_layout$position[2] > 0.5),
  "dynamic legend selector should avoid dense top-left corner"
)

external_legend_right <- choose_external_legend_layout(
  legend_nrow = 6,
  plot_width = 6.6,
  plot_height = 4.9
)
assert_true(
  identical(external_legend_right$position, "right"),
  "standard locus plot should place legend on the right outside the panel"
)

external_legend_bottom <- choose_external_legend_layout(
  legend_nrow = 7,
  plot_width = 5.2,
  plot_height = 4.4
)
assert_true(
  identical(external_legend_bottom$position, "bottom"),
  "narrow plot should place legend at the bottom outside the panel"
)

assign("lead_SNP", example_data$lead_snp, envir = .GlobalEnv)
assign("geneSymbol", example_data$phenotype_info, envir = .GlobalEnv)
assign("plink_bfile", NULL, envir = .GlobalEnv)
assign("trait", "SMOKE_GWAS", envir = .GlobalEnv)

plot_obj <- plot_qtl_association(
  qtl_all_chrom = "CHR.qtl",
  qtl_all_pvalue = "P.qtl",
  leadSNP_DF = example_data$merged_data,
  ld_df = NULL,
  gtf_path = gtf_path,
  region_recomb = NULL,
  recomb_path = recomb_path,
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
