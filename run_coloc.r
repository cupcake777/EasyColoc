#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# run_coloc.R 
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(yaml)
  library(parallel)
  library(tools)
  library(glue)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(jsonlite)
})

# ==============================================================================
# 1. Initialization
# ==============================================================================
message("[INIT] Loading Utility Modules...")
utils_files <- list.files("src", pattern = "^utils_.*\\.R$", full.names = TRUE)
sapply(utils_files, source)

cfg_global <- read_yaml("config/global.yaml")
cfg_gwas   <- read_yaml("config/gwas.yaml")
cfg_qtl    <- read_yaml("config/qtl.yaml")

# Read PP4 threshold from config (default to 0.8 if missing)
coloc_pp4_thresh <- if(!is.null(cfg_global$sig_threshold)) as.numeric(cfg_global$sig_threshold) else 0.8
message(glue("[INIT] Coloc PP4 Threshold set to: {coloc_pp4_thresh}"))

# Output Dirs
base_out_dir <- normalizePath(cfg_global$output_dir, mustWork = FALSE)
dir_abf <- file.path(base_out_dir, "abf")
dir_susie <- file.path(base_out_dir, "susie")
dir_plots <- file.path(base_out_dir, "plots")
for(d in c(base_out_dir, dir_abf, dir_susie, dir_plots)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

message(glue("[INIT] Output directories ready: {base_out_dir}"))

# Hash Table Init
hash_available <- FALSE
if (!is.null(cfg_global$hash_table_dir) && dir.exists(cfg_global$hash_table_dir)) {
    hash_available <- initialize_hash_system(cfg_global$hash_table_dir)
}

# ==============================================================================
# 2. Helpers (Locus Identification)
# ==============================================================================
identify_loci <- function(sumstats_dt, p_col, snp_col, chrom_col, pos_col,
                          plink_bfile = NULL, plink_bin = "plink",
                          sig_thresh = 1e-6, clump_kb = 1000, clump_r2 = 0.2) {

    message(glue("[LOCUS] Identifying loci (P < {sig_thresh})..."))
    if (!all(c(p_col, snp_col, chrom_col, pos_col) %in% names(sumstats_dt))) return(NULL)

    sig_dt <- sumstats_dt[sumstats_dt[[p_col]] <= sig_thresh, ]
    if (nrow(sig_dt) == 0) { warning("No significant SNPs."); return(NULL) }

    # PLINK Clumping (Preferred)
    if (!is.null(plink_bfile)) {
        message(glue("        Using PLINK clumping with: {basename(plink_bfile)}"))
        temp_assoc <- tempfile(fileext = ".qassoc"); temp_prefix <- tempfile()

        # Ensure unique SNPs for clumping input
        clump_input <- sig_dt[!duplicated(sig_dt[[snp_col]]), ]
        fwrite(clump_input[, .(SNP = get(snp_col), P = get(p_col))], temp_assoc, sep = "\t", quote = FALSE)

        cmd <- glue("{plink_bin} --bfile {plink_bfile} --clump {temp_assoc} ",
                    "--clump-p1 {sig_thresh} --clump-p2 {sig_thresh} ",
                    "--clump-r2 {clump_r2} --clump-kb {clump_kb} --out {temp_prefix}")

        system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
        clump_file <- paste0(temp_prefix, ".clumped")

        if (file.exists(clump_file)) {
            clump_res <- fread(clump_file)
            if (nrow(clump_res) > 0) {
                message(glue("[LOCUS] Found {nrow(clump_res)} clumps."))
                loci <- list()
                for(i in 1:nrow(clump_res)) {
                    lead_snp <- clump_res$SNP[i]
                    row <- sig_dt[sig_dt[[snp_col]] == lead_snp, ][1] # Take first if dup
                    if(nrow(row) > 0) {
                        loci[[i]] <- list(
                            chrom = as.character(row[[chrom_col]]),
                            pos = as.numeric(row[[pos_col]]),
                            snp = lead_snp,
                            p = as.numeric(row[[p_col]])
                        )
                    }
                }
                unlink(c(temp_assoc, paste0(temp_prefix, "*"))); return(loci)
            }
        }
        unlink(c(temp_assoc, paste0(temp_prefix, "*")))
    }

    # Fallback: Distance Pruning
    message("        Fallback to distance pruning...")
    sig_dt <- sig_dt[order(sig_dt[[p_col]]), ]
    loci <- list()
    while(nrow(sig_dt) > 0) {
        lead <- sig_dt[1, ]
        loci[[length(loci)+1]] <- list(chrom = as.character(lead[[chrom_col]]), pos = as.numeric(lead[[pos_col]]), snp = lead[[snp_col]], p = as.numeric(lead[[p_col]]))
        sig_dt <- sig_dt[!(sig_dt[[chrom_col]] == lead[[chrom_col]] & abs(sig_dt[[pos_col]] - lead[[pos_col]]) <= (clump_kb * 1000 / 2)), ]
    }
    return(loci)
}

# ==============================================================================
# 3. Main Execution
# ==============================================================================

run_pipeline <- function() {

  qtl_meta <- fread(cfg_qtl$qtl_info$file)
  message(glue("[DATA] Loaded {nrow(qtl_meta)} QTL datasets."))

  for (gwas_cfg in cfg_gwas$datasets) {
      separator <- strrep("=", 50)
      message(glue("\n{separator}\nProcessing GWAS: {gwas_cfg$name}\n{separator}"))

      if (!file.exists(gwas_cfg$file)) { warning(glue("File missing: {gwas_cfg$file}")); next }

      # ---------------------------------------------------------
      # 1. Load & Harmonize & LiftOver (to hg38)
      # ---------------------------------------------------------
      gwas_raw <- fread(gwas_cfg$file)
      gwas_std <- format_sumstats(gwas_raw, type="gwas", col_map=as.list(gwas_cfg$columns), case_control=(gwas_cfg$type=="cc"))

      # Harmonization handles LiftOver if source_build != target_build
      # We force target_build = "38" so all subsequent analysis is in hg38
      ref_fasta <- if(gwas_cfg$build=="hg19") cfg_global$ref_genome_hg19 else cfg_global$ref_genome_hg38

      gwas_harm <- run_gwaslab_harmonization(
           gwas_std,
           ref_fasta = ref_fasta,
           source_build = if(is.null(gwas_cfg$build)) "19" else gsub("hg", "", gwas_cfg$build),
           target_build = "38", 
           save_dir = cfg_global$harmonize_dir,
           dataset_id = gwas_cfg$id
      )
      
      # Ensure rsID column availability
      rsid_col_gwas <- NULL
      if ("rsid" %in% names(gwas_harm)) rsid_col_gwas <- "rsid"
      else if ("SNPID" %in% names(gwas_harm)) rsid_col_gwas <- "SNPID"

      plink_ref_hg38 <- cfg_global$plink_hg38

      loci_list <- identify_loci(
          gwas_harm,
          p_col="P", snp_col="SNPID", chrom_col="CHR", pos_col="POS",
          plink_bfile=plink_ref_hg38
          # sig_thresh = 1e-6 # Default used, change if needed
      )

      if (is.null(loci_list) || length(loci_list) == 0) { message("No significant loci."); next }

      # ---------------------------------------------------------
      # 3. Process Each Locus
      # ---------------------------------------------------------
      for (locus in loci_list) {
          message(glue("\n>>> Locus: {locus$snp} (hg38 chr{locus$chrom}:{locus$pos})"))

          # Locus is already hg38
          query_chrom <- locus$chrom
          query_pos <- locus$pos
          
          # Define Window (Flanking distance)
          # Use 500kb flank (total 1Mb window)
          flank_bp <- 500000 
          qtl_start <- max(1, query_pos - flank_bp)
          qtl_end   <- query_pos + flank_bp

          # Extract GWAS Subset for this region
          gwas_locus <- gwas_harm %>%
              filter(CHR == query_chrom & POS >= qtl_start & POS <= qtl_end)

          if(nrow(gwas_locus) == 0) next

          # Load Recomb Map (hg38) - Using utils_helpers.R function
          recomb_data <- load_recomb_map(query_chrom, qtl_start, qtl_end, cfg_global$recom)

          # Process QTLs
          process_qtl_wrapper <- function(row_idx) {
              tryCatch({
                  meta_row <- qtl_meta[row_idx, ]
                  qtl_id <- as.character(meta_row[[ cfg_qtl$qtl_info$columns$id ]])
                  qtl_file <- file.path(cfg_qtl$qtl_info$base_dir, meta_row[[ cfg_qtl$qtl_info$columns$all_filename ]])
                  qtl_n <- as.numeric(meta_row[[ cfg_qtl$qtl_info$columns$sample_size ]])

                  # Tabix Extract (hg38)
                  qtl_raw <- query_tabix_region(qtl_file, query_chrom, qtl_start, qtl_end)
                  if (is.null(qtl_raw) || nrow(qtl_raw) == 0) return(NULL)
                  if (ncol(qtl_raw) == length(cfg_qtl$QTL_all_header)) colnames(qtl_raw) <- cfg_qtl$QTL_all_header

                  pheno_col <- cfg_qtl$QTL_cols$phenotype
                  unique_phenos <- if (is.null(pheno_col)) "Combined" else unique(qtl_raw[[pheno_col]])

                  res_list <- list()

                  for(pheno in unique_phenos) {
                      qtl_sub <- if(pheno=="Combined") qtl_raw else qtl_raw[qtl_raw[[pheno_col]] == pheno, ]
                      qtl_std <- format_sumstats(qtl_sub, type="qtl", col_map=as.list(cfg_qtl$QTL_cols))

                      # Prepare QTL rsID
                      if(!is.null(rsid_col_gwas)) gwas_locus$rsid <- gwas_locus[[rsid_col_gwas]]
                      
                      # Try to find rsID in QTL
                      if ("variant_id" %in% names(qtl_std)) qtl_std$rsid <- qtl_std$variant_id # Fallback
                      if ("snp" %in% names(qtl_std) && any(grepl("^rs", qtl_std$snp[1:min(5,nrow(qtl_std))]))) qtl_std$rsid <- qtl_std$snp
                      if ("SNPID" %in% names(qtl_std) && any(grepl("^rs", qtl_std$SNPID[1:min(5,nrow(qtl_std))]))) qtl_std$rsid <- qtl_std$SNPID

                      # Merge
                      input_data <- prep_coloc_input_file(
                          gwas_locus, qtl_std,
                          list(pval="P", beta="BETA", se="SE", n="N"),
                          list(pval="P", beta="BETA", se="SE"),
                          use_hash_table = hash_available
                      )

                      if(nrow(input_data) < 30) next

                      # Analysis
                      res <- get_coloc_results(input_data, gwas_cfg$type, gwas_cfg$prop, NULL, qtl_n, 
                                               plink_bfile = plink_ref_hg38, plink_bin = "plink", use_susie=TRUE)

                      pp4 <- res$summary["PP.H4.abf"]

                      # Symbol Mapping
                      gene_sym <- pheno
                      try({
                          base <- gsub("\\..*", "", pheno)
                          s <- suppressMessages(bitr(base, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db"))
                          if(nrow(s) > 0) gene_sym <- s$SYMBOL[1]
                      }, silent=TRUE)

                      # Save SuSiE
                      if (!is.null(res$susie_result) && !is.null(res$susie_result$summary) && nrow(res$susie_result$summary)>0) {
                          fname <- gsub("[:/]", "_", glue("{gwas_cfg$id}_{qtl_id}_{gene_sym}_{locus$snp}_susie.csv"))
                          fwrite(as.data.frame(res$susie_result$summary), file.path(dir_susie, fname))
                      }

                      # Plot
                      if(!is.na(pp4) && pp4 > coloc_pp4_thresh) {
                          # Set globals for plot function
                          assign("lead_SNP", locus$snp, envir=.GlobalEnv)
                          assign("geneSymbol", gene_sym, envir=.GlobalEnv)
                          assign("tissue", qtl_id, envir=.GlobalEnv)
                          assign("plink_bfile", plink_ref_hg38, envir=.GlobalEnv)

                          tryCatch({
                              p <- plot_qtl_association(
                                  qtl_all_chrom="CHR.qtl", qtl_all_pvalue="P.qtl",
                                  leadSNP_DF=input_data,
                                  gtf_path=cfg_global$gene_anno,
                                  region_recomb=recomb_data
                              )
                              pname <- gsub("[:/]", "_", glue("{gwas_cfg$id}_{qtl_id}_{gene_sym}_{locus$snp}_coloc.pdf"))
                              ggsave(file.path(dir_plots, pname), p, width=10, height=6)
                          }, error=function(e) message(glue("Plot Error: {e$message}")))
                      }

                      res_list[[length(res_list)+1]] <- data.table(
                          Locus=locus$snp, Gene=gene_sym, PP4=pp4, n_snps=res$summary["nsnps"]
                      )
                  }
                  return(rbindlist(res_list))
              }, error=function(e) return(NULL))
          }

          locus_results <- mclapply(1:nrow(qtl_meta), process_qtl_wrapper, mc.cores = if(cfg_global$use_parallel) parallel::detectCores()-1 else 1)
          final_dt <- rbindlist(locus_results, fill=TRUE)

          if(nrow(final_dt) > 0) {
              final_dt <- final_dt[order(-PP4)]
              fwrite(final_dt, file.path(dir_abf, glue("{gwas_cfg$id}_{locus$snp}_results.csv")))
              message(glue("    -> Saved results for {locus$snp}"))
          }
      }
  }
}

run_pipeline()
