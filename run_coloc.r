#!/usr/bin/env Rscript
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
message("==============================================================")
message("EasyColoc v1.1 - Colocalization Analysis Pipeline")
message("==============================================================")

message("[INIT] Loading Utility Modules...")
utils_files <- list.files("src", pattern = "^utils_.*\\.R$", full.names = TRUE)
sapply(utils_files, source)

cfg_global <- read_yaml("config/global.yaml")
cfg_gwas   <- read_yaml("config/gwas.yaml")
cfg_qtl    <- read_yaml("config/qtl.yaml")

# Read colocalization analysis settings
coloc_pp4_thresh <- if(!is.null(cfg_global$coloc_settings$pp4_threshold)) {
    as.numeric(cfg_global$coloc_settings$pp4_threshold)
} else if(!is.null(cfg_global$sig_threshold)) {
    as.numeric(cfg_global$sig_threshold)
} else {
    0.8
}
susie_thresh <- if(!is.null(cfg_global$coloc_settings$susie_threshold)) {
    as.numeric(cfg_global$coloc_settings$susie_threshold)
} else {
    0.75
}
min_snps <- if(!is.null(cfg_global$coloc_settings$min_snps)) {
    as.integer(cfg_global$coloc_settings$min_snps)
} else {
    30
}
top_candidates <- if(!is.null(cfg_global$coloc_settings$top_candidates)) {
    as.integer(cfg_global$coloc_settings$top_candidates)
} else {
    100
}
maf_default <- if(!is.null(cfg_global$coloc_settings$maf_default)) {
    as.numeric(cfg_global$coloc_settings$maf_default)
} else {
    0.1
}
maf_na_replacement <- if(!is.null(cfg_global$coloc_settings$maf_na_replacement)) {
    as.numeric(cfg_global$coloc_settings$maf_na_replacement)
} else {
    0.05
}
maf_epsilon <- if(!is.null(cfg_global$coloc_settings$maf_epsilon)) {
    as.numeric(cfg_global$coloc_settings$maf_epsilon)
} else {
    1.0e-6
}
pvalue_floor <- if(!is.null(cfg_global$coloc_settings$pvalue_floor)) {
    as.numeric(cfg_global$coloc_settings$pvalue_floor)
} else {
    1.0e-300
}

# Read coloc priors from config
coloc_p1 <- if(!is.null(cfg_global$coloc_settings$p1)) as.numeric(cfg_global$coloc_settings$p1) else 1e-4
coloc_p2 <- if(!is.null(cfg_global$coloc_settings$p2)) as.numeric(cfg_global$coloc_settings$p2) else 1e-4
coloc_p12 <- if(!is.null(cfg_global$coloc_settings$p12)) as.numeric(cfg_global$coloc_settings$p12) else 5e-6

# Read harmonization settings

# Read harmonization settings
gwaslab_env <- if(!is.null(cfg_global$harmonization_settings$env_name)) {
    cfg_global$harmonization_settings$env_name
} else {
    "gwaslab"
}

# Read plot settings
plot_sig_threshold <- if(!is.null(cfg_global$plot_settings$significance_threshold)) {
    as.numeric(cfg_global$plot_settings$significance_threshold)
} else {
    5.0e-8
}
plot_window_bp <- if(!is.null(cfg_global$plot_settings$plot_window_bp)) {
    as.integer(cfg_global$plot_settings$plot_window_bp)
} else {
    200000
}
r2_breaks <- if(!is.null(cfg_global$plot_settings$r2_breaks)) {
    as.numeric(unlist(cfg_global$plot_settings$r2_breaks))
} else {
    c(0.2, 0.4, 0.6, 0.8, 1.0)
}
r2_colors <- if(!is.null(cfg_global$plot_settings$r2_colors)) {
    as.character(unlist(cfg_global$plot_settings$r2_colors))
} else {
    c("#313695", "#4575B4", "#74ADD1", "#FDB863", "#D73027")
}
lead_snp_color <- if(!is.null(cfg_global$plot_settings$lead_snp_color)) {
    cfg_global$plot_settings$lead_snp_color
} else {
    "#7F3C8D"
}

message(glue("[INIT] PP4 Threshold: {coloc_pp4_thresh}"))
message(glue("[INIT] SuSiE Threshold: {susie_thresh}"))
message(glue("[INIT] Min SNPs: {min_snps}"))
message(glue("[INIT] GWASLab Environment: {gwaslab_env}"))

base_out_dir <- normalizePath(cfg_global$output_dir, mustWork = FALSE)
dir_abf <- file.path(base_out_dir, "abf")
dir_susie <- file.path(base_out_dir, "susie")
dir_plots <- file.path(base_out_dir, "plots")
for(d in c(base_out_dir, dir_abf, dir_susie, dir_plots)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

message(glue("[INIT] Output: {base_out_dir}"))
hash_table <- TRUE
if (!is.null(cfg_global$hash_table_dir) && dir.exists(cfg_global$hash_table_dir)) {
    hash_table <- initialize_hash_system(cfg_global$hash_table_dir)
}

identify_loci <- function(sumstats_dt, p_col, snp_col, chrom_col, pos_col, plink_bfile = NULL, plink_bin = "plink",
                           clump_p1 = 5e-8, clump_p2 = 5e-8, clump_kb = 1000, clump_r2 = 0.1, dataset_id = NULL, keep_file = NULL) {
    message(glue("[LOCUS] Identifying loci (P < {clump_p1})..."))
    if (!all(c(p_col, snp_col, chrom_col, pos_col) %in% names(sumstats_dt))) {
        stop("Columns not found in sumstats_dt")
    }
    sig_dt <- sumstats_dt[sumstats_dt[[p_col]] <= clump_p1, ]
    if (nrow(sig_dt) == 0) { warning("No significant SNPs."); return(NULL) }
     if (!is.null(plink_bfile)) {
        message(glue("        Using PLINK clumping with: {basename(plink_bfile)}"))
        # DEBUG: Use non-temp files for debugging
        debug_dir <- "debug_plink"
        if (!dir.exists(debug_dir)) dir.create(debug_dir, recursive = TRUE)
        ds_id <- if(is.null(dataset_id)) "unknown" else dataset_id
        temp_assoc <- file.path(debug_dir, glue("clump_input_{ds_id}.qassoc"))
        temp_prefix <- file.path(debug_dir, glue("clump_output_{ds_id}"))

        clump_input <- sig_dt[!duplicated(sig_dt[[snp_col]]), ]
        clump_write <- clump_input[, c(snp_col, p_col), with=FALSE]
        colnames(clump_write) <- c("SNP", "P")
        fwrite(clump_write, temp_assoc, sep = "\t", quote = FALSE, col.names = TRUE)
        message(glue("        DEBUG: Wrote {nrow(clump_write)} SNPs to {temp_assoc}"))

        # Show first few lines of input
        message(glue("        DEBUG: Input sample: "))
        message(paste(readLines(temp_assoc)[1:5], collapse="\n        "))

        keep_cmd <- if (!is.null(keep_file) && file.exists(keep_file)) {
            glue("--keep {keep_file}")
        } else {
            ""
        }
        cmd <- glue("{plink_bin} --bfile {plink_bfile} --clump {temp_assoc} {keep_cmd} ",
                    "--clump-p1 {clump_p1} --clump-p2 {clump_p2} ",
                    "--clump-r2 {clump_r2} --clump-kb {clump_kb} --out {temp_prefix}")

        message(glue("        DEBUG: PLINK command: {cmd}"))
        sys_result <- system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)

        # Check log file
        log_file <- paste0(temp_prefix, ".log")
        if (file.exists(log_file)) {
            log_content <- readLines(log_file)
            message(glue("        DEBUG: PLINK log:"))
            message(paste(tail(log_content, 20), collapse="\n        "))
        }

        if (sys_result != 0) {
            warning(glue("[LOCUS] PLINK execution failed (exit code: {sys_result}). Check debug_plink/ directory."))
        } else {
            clump_file <- paste0(temp_prefix, ".clumped")
            if (file.exists(clump_file)) {
                clump_res <- fread(clump_file)
                if (nrow(clump_res) == 0) {
                    warning("[LOCUS] PLINK ran successfully but found 0 clumps. Possible SNP ID mismatch between sumstats and reference panel? Falling back to distance pruning.")
                } else {
                    message(glue("[LOCUS] Found {nrow(clump_res)} clumps via PLINK."))
                    loci <- list()
                    for(i in 1:nrow(clump_res)) {
                        lead_snp <- clump_res$SNP[i]
                        row <- sig_dt[sig_dt[[snp_col]] == lead_snp, ][1]
                        if(nrow(row) > 0) {
                            loci[[i]] <- list(
                                chrom = as.character(row[[chrom_col]]),
                                pos = as.numeric(row[[pos_col]]),
                                snp = lead_snp,
                                p = as.numeric(row[[p_col]])
                            )
                        }
                    }
                    return(loci)
                }
            }
        }
    }
    message("        Fallback to distance pruning...")
    sig_dt <- sig_dt[order(sig_dt[[p_col]]), ]
    loci <- list()
    while(nrow(sig_dt) > 0) {
        lead <- sig_dt[1, ]
        loci[[length(loci)+1]] <- list(
            chrom = as.character(lead[[chrom_col]]), 
            pos = as.numeric(lead[[pos_col]]), 
            snp = lead[[snp_col]], 
            p = as.numeric(lead[[p_col]])
        )
        lead_chrom <- as.character(lead[[chrom_col]])
        lead_pos <- as.numeric(lead[[pos_col]])
        radius_bp <- clump_kb * 1000
        sig_dt <- sig_dt[!(as.character(sig_dt[[chrom_col]]) == lead_chrom & 
                           abs(sig_dt[[pos_col]] - lead_pos) <= radius_bp), ]
    }
    return(loci)
}


run_pipeline <- function() {
  qtl_meta <- fread(cfg_qtl$qtl_info$file)
  message(glue("[DATA] Loaded {nrow(qtl_meta)} QTL datasets."))
  
  # Parallel logic: use n_cores if set, else detectCores()-1 if use_parallel=TRUE, else 1
  if (!is.null(cfg_global$n_cores)) {
      n_cores <- as.integer(cfg_global$n_cores)
  } else if (isTRUE(cfg_global$use_parallel)) {
      n_cores <- parallel::detectCores() - 1
  } else {
      n_cores <- 1
  }
  message(glue("[INIT] Using {n_cores} CPU cores"))
   for (gwas_cfg in cfg_gwas$datasets) {
       separator <- strrep("=", 50)
       message(glue("\n{separator}\nProcessing GWAS Dataset: {gwas_cfg$name}\n{separator}"))
        if (!file.exists(gwas_cfg$file)) { warning(glue("File missing: {gwas_cfg$file}")); next }
         gwas_raw <- fread(gwas_cfg$file)
          gwas_std <- format_sumstats(gwas_raw, type="gwas", col_map=as.list(gwas_cfg$columns),
                                   case_control=(gwas_cfg$type=="cc"))
          ref_fasta <- if(gwas_cfg$build=="hg19") cfg_global$ref_genome_hg19 else cfg_global$ref_genome_hg38
          # Get population first (needed for AF file naming)
          pop <- if(is.null(gwas_cfg$pop)) "EUR" else gwas_cfg$pop
          # Use pre-processed 1KG AF TSV for strand inference (no need for per-chromosome VCF)
          ref_vcf_1kg <- file.path(cfg_global[["1kg_af"]], glue("1KG_hg{gsub('hg', '', gwas_cfg$build)}_{pop}_AF.tsv.gz"))
          # Select dbSNP VCF for rsID annotation
          ref_dbsnp <- if(gwas_cfg$build=="hg19") cfg_global[["dbsnp_hg19"]] else cfg_global[["dbsnp_hg38"]]
          # Ref AF field is now统一为 "AF" in pre-processed TSV files
          ref_alt_freq <- "AF"
          # QTL build version from qtl.yaml
          qtl_build <- if(is.null(cfg_qtl$qtl_info$build)) "hg38" else cfg_qtl$qtl_info$build
            gwas_harm <- run_gwaslab_harmonization(
                gwas_std,
                ref_fasta = ref_fasta,
                ref_vcf = ref_vcf_1kg,
                ref_dbsnp = ref_dbsnp,
                ref_alt_freq = ref_alt_freq,
                source_build = if(is.null(gwas_cfg$build)) "19" else gsub("hg", "", gwas_cfg$build),
                target_build = gsub("hg", "", qtl_build),
                n_threads = cfg_global$n_cores,
                save_dir = cfg_global$harmonize_dir,
                dataset_id = gwas_cfg$id,
                env_name = gwaslab_env)
        if (!"N" %in% names(gwas_harm) && !is.null(gwas_cfg$sample_size_n)) {
            gwas_harm[, N := as.numeric(gwas_cfg$sample_size_n)]
        }
        # gwaslab的输出format是固定的
        plink_ref_hg38 <- cfg_global$plink_hg38
        rsid_col_gwas <- if("rsID" %in% names(gwas_harm)) "rsID" else if("SNPID" %in% names(gwas_harm)) "SNPID"
        p_col_locus <- if("P" %in% names(gwas_harm)) "P"
        chr_col_locus <- if("CHR" %in% names(gwas_harm)) "CHR"
        pos_col_locus <- if("POS" %in% names(gwas_harm)) "POS"
        if (is.null(p_col_locus) || is.null(rsid_col_gwas) || is.null(chr_col_locus) || is.null(pos_col_locus)) {
            stop("Missing required columns for locus identification")
        }
         # Get clump parameters: Dataset Config > Global Config > Hardcoded Defaults
         # p1 (primary threshold)
         clump_p1 <- if (!is.null(gwas_cfg$clump_p1)) {
             as.numeric(gwas_cfg$clump_p1)
         } else if (!is.null(gwas_cfg$clump_thres)) {
             as.numeric(gwas_cfg$clump_thres)  # backward compatibility
         } else if (!is.null(cfg_global$clump$p1)) {
             as.numeric(cfg_global$clump$p1)
         } else {
             5e-8  # hardcoded safety default
         }

         # p2 (secondary threshold) - defaults to p1 if not specified
         clump_p2 <- if (!is.null(gwas_cfg$clump_p2)) {
             as.numeric(gwas_cfg$clump_p2)
         } else if (!is.null(cfg_global$clump$p2)) {
             as.numeric(cfg_global$clump$p2)
         } else {
             clump_p1  # default to p1
         }

         # kb (window size)
         clump_kb <- if (!is.null(gwas_cfg$clump_kb)) {
             as.numeric(gwas_cfg$clump_kb)
         } else if (!is.null(cfg_global$clump$kb)) {
             as.numeric(cfg_global$clump$kb)
         } else {
             1000  # hardcoded safety default
         }

         # r2 (LD threshold)
         clump_r2 <- if (!is.null(gwas_cfg$clump_r2)) {
             as.numeric(gwas_cfg$clump_r2)
         } else if (!is.null(cfg_global$clump$r2)) {
             as.numeric(cfg_global$clump$r2)
         } else {
             0.1  # hardcoded safety default
         }

         message(glue("[LOCUS] Clump params: p1={clump_p1}, p2={clump_p2}, kb={clump_kb}, r2={clump_r2}"))
         # Determine keep_file based on population
         pop <- if(is.null(gwas_cfg$pop)) "EUR" else gwas_cfg$pop
         keep_base_dir <- dirname(cfg_global$plink_keep)
         keep_file <- file.path(keep_base_dir, glue("{pop}.sample"))
         if (!file.exists(keep_file)) {
             message(glue("        WARNING: Keep file not found: {keep_file}, using all samples"))
             keep_file <- NULL
         } else {
             message(glue("        Using population-specific LD: {pop}"))
         }

          loci_list <- identify_loci(
              gwas_harm,
              p_col = p_col_locus,
              snp_col = rsid_col_gwas,
              chrom_col = chr_col_locus,
              pos_col = pos_col_locus,
              plink_bfile = plink_ref_hg38,
              clump_p1 = clump_p1,
              clump_p2 = clump_p2,
              clump_kb = clump_kb,
              clump_r2 = clump_r2,
              dataset_id = gwas_cfg$id,
              keep_file = keep_file
          )
      if (is.null(loci_list) || length(loci_list) == 0) { message("No significant loci."); next }
       for (locus in loci_list) {
       message(glue("\n>>> Locus: {locus$snp} (hg38 chr{locus$chrom}:{locus$pos})"))
       query_chrom <- locus$chrom
       query_pos <- locus$pos
       flank_bp <- 500000
       qtl_start <- max(1, query_pos - flank_bp)
       qtl_end   <- query_pos + flank_bp
       gwas_locus <- gwas_harm %>%
       filter(CHR == query_chrom & POS >= qtl_start & POS <= qtl_end)
       if(nrow(gwas_locus) == 0) next
       recomb_data <- load_recomb_map(query_chrom, qtl_start, qtl_end, cfg_global$recom)
       process_qtl_wrapper <- function(row_idx) {
              tryCatch({
                  meta_row <- qtl_meta[row_idx, ]
                  qtl_id <- as.character(meta_row[[ cfg_qtl$qtl_info$columns$id ]])
                  qtl_file <- meta_row[[ cfg_qtl$qtl_info$columns$all_filename ]]
                  qtl_n <- as.numeric(meta_row[[ cfg_qtl$qtl_info$columns$sample_size ]])
                  qtl_raw <- query_tabix_region(qtl_file, query_chrom, qtl_start, qtl_end)
                  if (is.null(qtl_raw) || nrow(qtl_raw) == 0) return(NULL)
                  if (ncol(qtl_raw) == length(cfg_qtl$QTL_all_header)) colnames(qtl_raw) <- cfg_qtl$QTL_all_header
                  pheno_col <- cfg_qtl$QTL_cols$phenotype
                  unique_phenos <- if (is.null(pheno_col)) "Combined" else unique(qtl_raw[[pheno_col]])
                  res_list <- list()
                  for(pheno in unique_phenos) {
                      qtl_sub <- if(pheno=="Combined") qtl_raw else qtl_raw[qtl_raw[[pheno_col]] == pheno, ]
                      qtl_std <- format_sumstats(qtl_sub, type="qtl", col_map=as.list(cfg_qtl$QTL_cols))
                      if(!is.null(rsid_col_gwas)) gwas_locus$rsid <- gwas_locus[[rsid_col_gwas]]
                      if ("variant_id" %in% names(qtl_std)) qtl_std$rsid <- qtl_std$variant_id
                        input_data <- prep_coloc_input_file(
                            gwas_locus, qtl_std,
                            list(pval="P.gwas", beta="BETA.gwas", se="SE.gwas", n="N.gwas"),
                            list(pval="P.qtl", beta="BETA.qtl", se="SE.qtl"),
                            use_hash_table = hash_table,
                            min_snps = min_snps,
                            pvalue_floor = pvalue_floor
                        )
                        if(nrow(input_data) < min_snps) next
                        res <- get_coloc_results(
                            input_data, gwas_cfg$type, gwas_cfg$prop, gwas_cfg$sample_size_n, qtl_n,
                            plink_bfile = plink_ref_hg38,
                            plink_bin = "plink",
                            use_susie = TRUE,
                            susie_threshold = susie_thresh,
                            maf_default = maf_default,
                            maf_na_replacement = maf_na_replacement,
                            maf_epsilon = maf_epsilon,
                            keep_file = cfg_global$plink_keep,
                            p1 = coloc_p1,
                            p2 = coloc_p2,
                            p12 = coloc_p12
                        )
                        pp4 <- as.numeric(res$summary["PP.H4.abf"])
                        gene_sym_for_filename <- pheno
                        try({
                            base <- gsub("\\..*", "", pheno)
                            s <- suppressMessages(bitr(base, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db"))
                            if(nrow(s) > 0) gene_sym_for_filename <- s$SYMBOL[1]
                        }, silent=TRUE)
                           if(!is.na(pp4) && pp4 > coloc_pp4_thresh) {
                               if (!is.null(res$susie_result) && !is.null(res$susie_result$summary) && nrow(res$susie_result$summary)>0) {
                                  fname <- sanitize_filename(glue("{gwas_cfg$id}_{qtl_id}_{gene_sym_for_filename}_susie.csv"))
                                  fwrite(as.data.frame(res$susie_result$summary), file.path(dir_susie, fname))
                               }
                           }

                          if(!is.na(pp4) && pp4 > coloc_pp4_thresh) {
                             tryCatch({
                                 assign("lead_SNP", locus$snp, envir = .GlobalEnv)
                                 assign("geneSymbol", gene_sym_for_filename, envir = .GlobalEnv)
                                 assign("plink_bfile", plink_ref_hg38, envir = .GlobalEnv)
                                 assign("trait", gwas_cfg$id, envir = .GlobalEnv)
                                  p <- plot_qtl_association(
                                      qtl_all_chrom = "CHR.qtl",
                                      qtl_all_pvalue = "P.qtl",
                                      leadSNP_DF = input_data,
                                      ld_df = NULL,
                                      gtf_path = cfg_global$gene_anno,
                                      region_recomb = recomb_data,
                                      show_lead_line = FALSE,
                                      plot_window_bp = plot_window_bp,
                                      significance_threshold = plot_sig_threshold,
                                      r2_breaks = r2_breaks,
                                      r2_colors = r2_colors,
                                      lead_snp_color = lead_snp_color
                                  )
                                  pname_pdf <- sanitize_filename(glue("{gwas_cfg$id}_{qtl_id}_{gene_sym_for_filename}_coloc.pdf"))
                                  pname_png <- sanitize_filename(glue("{gwas_cfg$id}_{qtl_id}_{gene_sym_for_filename}_coloc.png"))
                                  pdf_path <- file.path(dir_plots, pname_pdf)
                                  png_path <- file.path(dir_plots, pname_png)
                                  pdf_saved <- tryCatch({
                                    ggsave(pdf_path, p, width=12, height=10, device = cairo_pdf)
                                    TRUE
                                  }, error = function(e) {
                                    FALSE
                                  })

                                  if (!pdf_saved) {
                                    ggsave(png_path, p, width=10, height=6, dpi=300)
                                    message(glue("[PLOT] Saved as PNG: {pname_png}"))
                                  } else {
                                    message(glue("[PLOT] Saved as PDF: {pname_pdf}"))
                                  }
                            }, error=function(e) message(glue("Plot Error: {e$message}")))
                       }

                         res_list[[length(res_list)+1]] <- data.table(
                             GWAS_ID=gwas_cfg$id, QTL_ID=qtl_id, Locus=locus$snp, Gene=gene_sym_for_filename, PP4=pp4, n_snps=res$summary["nsnps"]
                         )
                  }
                  return(rbindlist(res_list))
              }, error=function(e) return(NULL))
          }
           locus_results <- mclapply(1:nrow(qtl_meta), process_qtl_wrapper, mc.cores = n_cores)
           final_dt <- rbindlist(locus_results, fill=TRUE)
             if(nrow(final_dt) > 0) {
                  final_dt <- final_dt[order(-PP4)]
                  results_fname <- sanitize_filename(glue("{gwas_cfg$id}_{locus$snp}_locus_results.csv"))
                  fwrite(final_dt, file.path(dir_abf, results_fname))
                  message(glue("-> Saved results for {locus$snp}"))
              }
       }
   }
   message("\n==============================================================")
   message("[SUM] Summary all results...")
   message("==============================================================")
   tryCatch({
       merge_all_results(
           output_dir = base_out_dir,
           pp4_threshold = coloc_pp4_thresh,
           merge_susie = TRUE,
           save_summary = TRUE
       )
   }, error = function(e) {
       warning(glue("[SUM] Failed to merge results: {e$message}"))
   })

}
run_pipeline()
message("==============================================================")
message("EasyColoc Analysis Complete!")
message("==============================================================")
