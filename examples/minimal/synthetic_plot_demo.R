#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

utils_files <- list.files("src", pattern = "^utils_.*\\.R$", full.names = TRUE)
invisible(lapply(utils_files, source))

cfg_global <- yaml::read_yaml("config/global.yaml")
gtf_path <- if (!is.null(cfg_global$gene_anno) && file.exists(cfg_global$gene_anno)) {
  cfg_global$gene_anno
} else {
  NULL
}
recomb_path <- if (!is.null(cfg_global$recom) && file.exists(cfg_global$recom)) {
  cfg_global$recom
} else {
  NULL
}

make_synthetic_locus <- function() {
  set.seed(20260327)

  n_bg <- 80
  bg_pos <- sort(sample(seq(43120000, 43580000, by = 500), n_bg))
  bg_p <- runif(n_bg, min = 5e-3, max = 0.08)

  signal_pos <- c(43562000, 43564000, 43566000, 43568000, 43570000, 43572000)
  signal_p <- c(2e-6, 8e-7, 2e-8, 7e-8, 3e-6, 9e-6)
  signal_rs <- c("rs2819340", "rs3791138", "rs10890255", "rs3791137", "rs2004899", "rs2819341")

  data.frame(
    rsid = c(paste0("rsBG", seq_len(n_bg)), signal_rs),
    CHR.qtl = "chr1",
    POS.qtl = c(bg_pos, signal_pos),
    P.qtl = c(bg_p, signal_p),
    P.gwas = c(runif(n_bg, min = 1e-3, max = 0.2), c(3e-5, 6e-6, 1e-9, 5e-6, 9e-5, 7e-5)),
    stringsAsFactors = FALSE
  )
}

synthetic_df <- make_synthetic_locus()
credible_set <- data.frame(
  snp = c("rs2819340", "rs3791138", "rs2004899"),
  SNP.PP.H4 = c(0.42, 0.31, 0.18),
  stringsAsFactors = FALSE
)

assign("lead_SNP", "rs10890255", envir = .GlobalEnv)
assign("geneSymbol", "ENST00000654683.1|CCDC30|chr1:42482663-42484158|+", envir = .GlobalEnv)
assign("plink_bfile", NULL, envir = .GlobalEnv)
assign("trait", "SMOKE_GWAS", envir = .GlobalEnv)

plot_obj <- plot_qtl_association(
  qtl_all_chrom = "CHR.qtl",
  qtl_all_pvalue = "P.qtl",
  leadSNP_DF = synthetic_df,
  ld_df = NULL,
  gtf_path = gtf_path,
  region_recomb = NULL,
  recomb_path = recomb_path,
  show_lead_line = FALSE,
  qtl_type = "postnatal",
  phenotype_info = "ENST00000654683.1|CCDC30|chr1:42482663-42484158|+",
  credible_set = credible_set,
  plot_window_bp = 200000
)

output_dir <- "examples/minimal/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save_res <- save_plot_with_fallback(
  plot_obj = plot_obj,
  pdf_path = file.path(output_dir, "synthetic_locus_demo.pdf"),
  png_path = file.path(output_dir, "synthetic_locus_demo.png"),
  plot_width = 6.6,
  plot_height = 4.9,
  png_scale = 1,
  png_dpi = 240
)

message("[DEMO] Saved minimal synthetic panel to: ", save_res$path)
