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

resolve_synthetic_plot_inputs <- function(cfg_global, fixture_dir = "examples") {
  list(
    gtf_path = resolve_existing_file(c(
      file.path(fixture_dir, "smoke_hg38_chr1.gtf"),
      cfg_global$gene_anno
    )),
    recomb_path = resolve_recomb_prefix(c(
      file.path(fixture_dir, "CHB"),
      cfg_global$recom
    ))
  )
}

make_synthetic_locus <- function() {
  set.seed(20260327)

  n_bg <- 80
  bg_pos <- sort(sample(seq(43120000, 43580000, by = 500), n_bg))
  bg_p <- runif(n_bg, min = 5e-3, max = 0.08)

  signal_pos <- c(43562000, 43564000, 43566000, 43568000, 43570000, 43572000)
  signal_p <- c(2e-6, 8e-7, 2e-8, 7e-8, 3e-6, 9e-6)
  signal_rs <- c("rs2819340", "rs3791138", "rs10890255", "rs3791137", "rs2004899", "rs2819341")

  df <- data.frame(
    rsid = c(paste0("rsBG", seq_len(n_bg)), signal_rs),
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
