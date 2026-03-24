# ------------------------------------------------------------------------------
# src/utils_format.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(glue)
  library(jsonlite)
  library(rlang)
})

`%notin%` <- Negate(`%in%`)
if (file.exists("src/utils_hash.R")) {
  source("src/utils_hash.R")
} else {
  warning("Hash conversion script not exist!.")
}
add_variant_id <- function(df, chr_col = "CHR", pos_col = "POS", ea_col = "EA", nea_col = "NEA") {
    req_cols <- c(chr_col, pos_col, ea_col, nea_col)
    if (!all(req_cols %in% names(df))) {
        stop(glue("Missing columns for ID creation: {paste(setdiff(req_cols, names(df)), collapse=', ')}"))
    }
    v_id <- paste0("chr", df[[chr_col]], "_", df[[pos_col]], "_", df[[ea_col]], "_", df[[nea_col]])
    v_id_flip <- paste0("chr", df[[chr_col]], "_", df[[pos_col]], "_", df[[nea_col]], "_", df[[ea_col]])
    if (inherits(df, "data.table")) {
        df[, variant_id := v_id]
        df[, variant_id_flip := v_id_flip]
    } else {
        df$variant_id <- v_id
        df$variant_id_flip <- v_id_flip
    }
    
    return(df)
}

format_sumstats <- function(df, type, col_map, case_control = FALSE) {
  df_std <- as.data.table(copy(df))

  target_names <- c(
    snp  = "SNPID",
    chrom = "CHR",
    pos  = "POS",
    a1   = "EA",
    a2   = "NEA",
    beta = "BETA",
    se   = "SE",
    pval = "P",
    n    = "N",
    maf  = "EAF",
    af   = "EAF"
  )

  for (key in names(col_map)) {
    if (key %in% names(target_names) && col_map[[key]] %in% names(df_std)) {
      setnames(df_std, old = col_map[[key]], new = target_names[[key]])
    }
  }

  numeric_cols <- c("BETA", "SE", "P", "POS", "N", "EAF")
  for (col in numeric_cols) {
    if (col %in% names(df_std)) df_std[[col]] <- as.numeric(df_std[[col]])
  }
  if (case_control && "BETA" %notin% names(df_std) && "OR" %in% names(df_std)) {
    df_std[, BETA := log(OR)]
  }

  return(df_std)
}

run_gwaslab_harmonization <- function(sumstats_dt,
                                      ref_fasta,
                                      ref_vcf = NULL,
                                      ref_dbsnp = NULL,
                                      ref_alt_freq = "AF",
                                      source_build = "19",
                                      target_build = "38",
                                      n_threads = 4,
                                      env_name = "gwaslab",
                                      save_dir = NULL,
                                      dataset_id = NULL,
                                      input_file = NULL,
                                      verbose = TRUE) {

  log_message <- function(..., verbose = TRUE) {
    if (isTRUE(verbose)) {
      message(...)
    }
  }
  log_message("--- Starting Harmonization & LiftOver ---", verbose = verbose)

  if (!file.exists(ref_fasta)) {
    warning(glue("Reference FASTA not found: {ref_fasta}. Skipping harmonization."))
    return(sumstats_dt)
  }

  source_file <- input_file

  use_cache <- !is.null(save_dir) && !is.null(dataset_id)
  final_output_file <- NULL

  if (use_cache) {
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
      ref_name <- tools::file_path_sans_ext(basename(ref_fasta))
      # Cache filename includes build info
      final_output_file <- file.path(save_dir, glue("{dataset_id}_b{source_build}to{target_build}_harmonized.tsv"))

      if (file.exists(final_output_file)) {
          cache_info <- file.info(final_output_file)
          cache_mtime <- cache_info$mtime
          cache_size_ok <- !is.na(cache_info$size) && cache_info$size > 0
          ref_mtime <- file.info(ref_fasta)$mtime
          input_mtime <- if (!is.null(source_file) && file.exists(source_file)) {
              file.info(source_file)$mtime
          } else {
              NA
          }
          dep_mtime <- max(c(ref_mtime, input_mtime), na.rm = TRUE)
          cache_ok <- cache_size_ok && !is.na(cache_mtime) && !is.na(dep_mtime) && cache_mtime >= dep_mtime
          if (cache_ok) {
              log_message(glue("Found cached harmonized file: {final_output_file}"), verbose = verbose)
              log_message("Loading from cache...", verbose = verbose)
              return(fread(final_output_file))
          }
          log_message(glue("Cached file is outdated or empty: {final_output_file}"), verbose = verbose)
      }
  }

  temp_id <- format(Sys.time(), "%H%M%S")
  input_tmp_file <- glue("temp_gwas_in_{temp_id}.tsv")
  output_file <- if(use_cache) final_output_file else glue("temp_gwas_out_{temp_id}.tsv")
  py_script_file <- glue("run_gwaslab_{temp_id}.py")

  on.exit({
    if (file.exists(input_tmp_file)) file.remove(input_tmp_file)
    if (file.exists(py_script_file)) file.remove(py_script_file)
    if (!use_cache && file.exists(output_file)) file.remove(output_file)
  }, add = TRUE)

   fwrite(sumstats_dt, input_tmp_file, sep = "\t", na = "NA", quote = FALSE)

    # Build VCF path pattern for strand inference
     vcf_pattern <- if (!is.null(ref_vcf)) ref_vcf else "NULL"
     dbsnp_pattern <- if (!is.null(ref_dbsnp)) ref_dbsnp else "NULL"

     # Build Python script as a single string to avoid paste0 quote issues
     py_script <- sprintf('
import gwaslab as gl
import pandas as pd
import numpy as np
import sys

try:
    # Loading sumstats with GWASLab...
    mysumstats = gl.Sumstats("%s", fmt="gwaslab", build="%s")
    mysumstats.basic_check()

    # Clean numeric columns to prevent type errors in GWASLab
    for col in mysumstats.data.columns:
        if mysumstats.data[col].dtype == "object":
            numeric_converted = pd.to_numeric(mysumstats.data[col], errors="coerce")
            if numeric_converted.notna().mean() > 0.1:
                mysumstats.data[col] = numeric_converted

    ref_seq = "%s"
    ref_rsid_vcf = "%s"
    ref_infer = "%s"  # Direct path to TSV file (no {chr} pattern needed)
    ref_af = "%s"
    src_b = "%s"
    tgt_b = "%s"

    # Use GWASLab v4 harmonize() function
    if ref_rsid_vcf != "NULL" or ref_infer != "NULL":
        # Running GWASLab harmonization...
        mysumstats.harmonize(
            ref_seq=ref_seq,
            ref_rsid_vcf=ref_rsid_vcf if ref_rsid_vcf != "NULL" else None,
            ref_infer=ref_infer,
            ref_alt_freq=ref_af,
            sweep_mode=True,
            threads=%d
        )

    # Perform LiftOver if needed
    if src_b != tgt_b:
        # Performing LiftOver...
        mysumstats.liftover(from_build=src_b, to_build=tgt_b, remove=True)

    mysumstats.data.to_csv("%s", sep="\\t", index=False, float_format="%%.6g")

except Exception as e:
    # Error in GWASLab
    import traceback
    traceback.print_exc()
    sys.exit(1)
 ', input_tmp_file, source_build,
   ref_fasta, dbsnp_pattern, vcf_pattern, ref_alt_freq, source_build, target_build,
   n_threads,
   output_file)

     writeLines(py_script, py_script_file)

  cmd <- glue("micromamba run -n {env_name} python {py_script_file}")
   log_message("Running Python script...", verbose = verbose)
  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

   if (exit_code == 0) {
       if (file.exists(output_file)) {
         res_dt <- fread(output_file)
          log_message(glue("Processing complete. Input SNPs: {nrow(sumstats_dt)} -> Output SNPs: {nrow(res_dt)}"), verbose = verbose)
          log_message(glue("  Output columns: {paste(names(res_dt), collapse=', ')}"), verbose = verbose)
          
          # Post-harmonization validation
          n_invalid <- 0
          if ("P" %in% names(res_dt)) {
              n_invalid_p <- sum(res_dt$P < 0 | res_dt$P > 1 | is.na(res_dt$P))
              if (n_invalid_p > 0) {
                  warning(glue("Found {n_invalid_p} invalid P-values (< 0, > 1, or NA) after harmonization"))
                  n_invalid <- n_invalid + n_invalid_p
              }
          }
          if ("SE" %in% names(res_dt)) {
              n_invalid_se <- sum(res_dt$SE <= 0 | is.na(res_dt$SE))
              if (n_invalid_se > 0) {
                  warning(glue("Found {n_invalid_se} invalid SE values (<= 0 or NA) after harmonization"))
                  n_invalid <- n_invalid + n_invalid_se
              }
          }
          if ("BETA" %in% names(res_dt) && "SE" %in% names(res_dt) && "P" %in% names(res_dt)) {
              # Check BETA/SE/P consistency (Z-score should roughly match)
              res_dt$Z_from_beta <- res_dt$BETA / res_dt$SE
              res_dt$Z_from_p <- qnorm(res_dt$P / 2, lower.tail = FALSE) * sign(res_dt$BETA)
              res_dt$Z_diff <- abs(res_dt$Z_from_beta - res_dt$Z_from_p)
              n_inconsistent <- sum(res_dt$Z_diff > 1.0, na.rm = TRUE)  # Allow ~1 Z-score unit tolerance
              if (n_inconsistent > nrow(res_dt) * 0.05) {  # >5% inconsistent = problem
                  warning(glue("Found {n_inconsistent}/{nrow(res_dt)} ({round(n_inconsistent/nrow(res_dt)*100, 1)}%) variants with inconsistent BETA/SE/P after harmonization"))
              }
              res_dt$Z_from_beta <- NULL
              res_dt$Z_from_p <- NULL
              res_dt$Z_diff <- NULL
          }
          if (n_invalid > 0) {
              warning(glue("Total post-harmonization QC failures: {n_invalid} issues detected"))
          } else {
              log_message("Post-harmonization validation: PASS", verbose = verbose)
          }
          if(use_cache) log_message(glue("Saved result to: {output_file}"), verbose = verbose)
         return(res_dt)
       } else {
         stop("GWASLab output file missing.")
       }
   } else {
     stop("GWASLab processing failed! Check logs above.")
   }
}

# format input
prep_coloc_input_file <- function(gwas_df, qtl_df,
                                  gwas_cols = list(pval="P", beta="BETA", se="SE", n="N"),
                                  qtl_cols = list(pval="pval_nominal", beta="slope", se="slope_se"),
                                  use_hash_table = TRUE,
                                  min_snps = 30,
                                  pvalue_floor = 1e-300,
                                  verbose = FALSE) {

    # localized message function to control verbosity
    msg <- function(...) {
        if(verbose) message(...)
    }

    msg("\n--- MERGING DATA ---")
    gwas_df <- as.data.table(gwas_df)
    qtl_df <- as.data.table(qtl_df)
    msg(glue("[MERGE] Input: GWAS={nrow(gwas_df)} SNPs, QTL={nrow(qtl_df)} SNPs"))
    # Standardize column names for merging
    # Ensure GWAS has 'rsid' (lowercase), QTL has 'variant_id' (lowercase)
    # First remove existing 'rsid' if 'rsID' exists, then rename
    if ("rsID" %in% names(gwas_df)) {
        if ("rsid" %in% names(gwas_df)) {
            gwas_df[, rsid := NULL]
        }
        setnames(gwas_df, "rsID", "rsid")
    }
    if ("variant_id" %in% names(qtl_df) && !"variant_id" %in% names(qtl_df)) {
        # Already lowercase
    }
    # Check for rsID columns
    gwas_has_rsid <- any(tolower(names(gwas_df)) %in% c("snpid", "snp", "rsid"))
    qtl_has_rsid <- any(tolower(names(qtl_df)) %in% c("snpid", "variant_id", "snp", "rsid"))
    qtl_has_pos <- all(c("CHR", "POS") %in% names(qtl_df)) ||
                   all(c("chr", "end") %in% names(qtl_df))
    gwas_has_pos <- all(c("CHR", "POS") %in% names(gwas_df))
    gwas_has_allele <- all(c("EA", "NEA") %in% names(gwas_df))
    qtl_has_allele <- all(c("ref", "alt") %in% names(qtl_df))
    if ("chr" %in% names(qtl_df) && "CHR" %notin% names(qtl_df)) {
        setnames(qtl_df, "chr", "CHR")
    }
    if ("end" %in% names(qtl_df) && "POS" %notin% names(qtl_df)) {
        setnames(qtl_df, "end", "POS")
    }
    if ("ref" %in% names(qtl_df) && "EA" %notin% names(qtl_df)) {
        if ("alt" %in% names(qtl_df)) {
            setnames(qtl_df, c("ref", "alt"), c("REF", "ALT"))
        } 
    }
    if ("REF" %in% names(qtl_df) && "EA" %notin% names(qtl_df)) {
         qtl_df[, `:=`(EA = ALT, NEA = REF)]
    }
    merged_dt <- NULL
    merge_strategy <- "none"
        msg("[MERGE] Attempting Strategy 1: rsID direct match...")
        # Use rsID format for matching (after standardization)
        gwas_id_col <- if ("rsid" %in% names(gwas_df)) "rsid" else if ("SNPID" %in% names(gwas_df)) "SNPID" else if ("SNP" %in% names(gwas_df)) "SNP" else "rsid"
        qtl_id_col <- if ("variant_id" %in% names(qtl_df)) "variant_id" else if ("SNPID" %in% names(qtl_df)) "SNPID" else if ("SNP" %in% names(qtl_df)) "SNP" else "snp"

        # Rename duplicate columns BEFORE merging
        dup_cols <- c("CHR", "POS", "EA", "NEA", "BETA", "SE", "P", "N", "EAF", "AF")
        for (col in dup_cols) {
            if (col %in% names(gwas_df) && col %in% names(qtl_df)) {
                new_name_qtl <- paste0(col, ".qtl")
                new_name_gwas <- paste0(col, ".gwas")
                if (!(new_name_qtl %in% names(qtl_df))) {
                    setnames(qtl_df, col, new_name_qtl)
                }
                if (!(new_name_gwas %in% names(gwas_df))) {
                    setnames(gwas_df, col, new_name_gwas)
                }
            }
        }

        # Manual merge to avoid data.table segfault
        common_snps <- intersect(gwas_df[[gwas_id_col]], qtl_df[[qtl_id_col]])
        if (length(common_snps) > 0) {
            gwas_subset <- gwas_df[gwas_df[[gwas_id_col]] %in% common_snps, ]
            qtl_subset <- qtl_df[qtl_df[[qtl_id_col]] %in% common_snps, ]
            merged_dt <- merge(
                gwas_subset,
                qtl_subset,
                by.x = gwas_id_col,
                by.y = qtl_id_col,
                all = FALSE
            )
            merged_dt <- as.data.table(merged_dt)
            merge_strategy <- "rsid_direct"
            # Priority: rsid.gwas > SNPID.gwas > SNP.gwas > SNPID > SNP
            # This is critical for SuSiE which needs rsIDs to query PLINK
            if ("rsid.gwas" %in% names(merged_dt)) {
                merged_dt$snp <- merged_dt$rsid.gwas
            } else if ("SNPID.gwas" %in% names(merged_dt)) {
                merged_dt$snp <- merged_dt$SNPID.gwas
            } else if ("SNP.gwas" %in% names(merged_dt)) {
                merged_dt$snp <- merged_dt$SNP.gwas
            } else if ("rsid" %in% names(merged_dt)) {
                merged_dt$snp <- merged_dt$rsid
            } else if ("SNPID" %in% names(merged_dt)) {
                merged_dt$snp <- merged_dt$SNPID
            } else if ("SNP" %in% names(merged_dt)) {
                merged_dt$snp <- merged_dt$SNP
            } else {
                merged_dt$snp <- merged_dt[[gwas_id_col]]
            }
        } else {
            msg("[MERGE] ✗ Strategy 1 failed: No rsID overlap")
            merged_dt <- NULL
    }
    if ((is.null(merged_dt) || nrow(merged_dt) < min_snps) && use_hash_table && exists("convert_rsid_to_pos")) {
        msg("[MERGE] Attempting Strategy 2: hash_table conversion...")
        gwas_id_col <- if ("rsid" %in% names(gwas_df)) "rsid" else if ("SNPID" %in% names(gwas_df)) "SNPID" else if ("SNP" %in% names(gwas_df)) "SNP" else NULL
        if (!is.null(gwas_id_col) && gwas_has_pos) {
            gwas_chrpos <- convert_rsid_to_pos(
                gwas_df[[gwas_id_col]],
                gwas_df$CHR
            )
            gwas_df$gwas_chrpos_id <- gwas_chrpos
            n_converted <- sum(!is.na(gwas_chrpos))
            msg(glue("[MERGE]   Converted {n_converted}/{nrow(gwas_df)} GWAS rsIDs ({round(n_converted/nrow(gwas_df)*100, 1)}%)"))
            qtl_id_col <- if ("variant_id" %in% names(qtl_df)) "variant_id" else if ("snp" %in% names(qtl_df)) "snp" else NULL
            if (!is.null(qtl_id_col)) {
                # Use dplyr to avoid segfault
                merged_dt <- inner_join(
                    gwas_df,
                    qtl_df,
                    by = c("gwas_chrpos_id" = qtl_id_col),
                    suffix = c(".gwas", ".qtl"),
                    relationship = "many-to-many"
                )
                merged_dt <- as.data.table(merged_dt)
                     msg(glue("[MERGE] ✓ Strategy 2: {nrow(merged_dt)} SNPs matched"))
                     merge_strategy <- "hash_conversion"
                     # Priority: rsid.gwas > SNPID.gwas > SNP.gwas > rsid > SNP > gwas_chrpos_id
                     rsid_col <- if ("rsid.gwas" %in% names(merged_dt)) "rsid.gwas"
                                 else if ("SNPID.gwas" %in% names(merged_dt)) "SNPID.gwas"
                                 else if ("SNP.gwas" %in% names(merged_dt)) "SNP.gwas"
                                 else if ("rsid" %in% names(merged_dt)) "rsid"
                                 else if ("SNPID" %in% names(merged_dt)) "SNPID"
                                 else if ("SNP" %in% names(merged_dt)) "SNP"
                                 else NULL
                     if (!is.null(rsid_col)) {
                         merged_dt$snp <- merged_dt[[rsid_col]]
                     } else {
                         setnames(merged_dt, "gwas_chrpos_id", "snp")
                     }
                     snp_sample <- head(merged_dt$snp[!is.na(merged_dt$snp)], 3)
                     if (length(snp_sample) > 0 && !all(grepl("^rs", snp_sample))) {
                         msg(glue("[MERGE]   Warning: Non-rsID format detected: {paste(snp_sample, collapse=', ')}"))
                     }
                 } else {
                     msg("[MERGE] ✗ Strategy 2 failed: No matches after hash conversion")
                 }
            }
        }
    if ((is.null(merged_dt) || nrow(merged_dt) < min_snps) && gwas_has_pos && qtl_has_pos) {
        msg("[MERGE] Attempting Strategy 3: position + allele matching...")
        if (gwas_has_allele) {
            gwas_df <- add_variant_id(gwas_df, "CHR", "POS", "EA", "NEA")
        }
        qtl_ref <- if("REF" %in% names(qtl_df)) "REF" else "EA"
        qtl_alt <- if("ALT" %in% names(qtl_df)) "ALT" else "NEA"
        if (qtl_ref %in% names(qtl_df) && qtl_alt %in% names(qtl_df)) {
             qtl_df <- add_variant_id(qtl_df, "CHR", "POS", qtl_ref, qtl_alt)
        }
        if ("variant_id" %in% names(gwas_df) && "variant_id" %in% names(qtl_df)) {
            # Forward Match
            merged_forward <- inner_join(
                gwas_df, qtl_df,
                by = "variant_id",
                suffix = c(".gwas", ".qtl"),
                relationship = "many-to-many"
            )
            merged_forward <- as.data.table(merged_forward)
            merged_flip <- inner_join(
                gwas_df, qtl_df,
                by = c("variant_id_flip" = "variant_id"),
                suffix = c(".gwas", ".qtl"),
                relationship = "many-to-many"
            )
            merged_flip <- as.data.table(merged_flip)
            # Mark flipped SNPs and flip BETA values
            # CRITICAL: Check if gwaslab already flipped these variants to prevent double-flipping
            if (nrow(merged_flip) > 0) {
                merged_flip$allele_flipped <- TRUE
                # Check if gwaslab ALLELE_FLIPPED column exists (indicates prior harmonization flips)
                if ("ALLELE_FLIPPED" %in% names(merged_flip)) {
                    # gwaslab tracked flips - check if any were actually flipped
                    already_flipped <- any(merged_flip$ALLELE_FLIPPED == TRUE, na.rm = TRUE)
                    if (already_flipped) {
                        warning("Detected ALLELE_FLIPPED=TRUE from gwaslab - variants already flipped during harmonization. Skipping downstream flip to prevent double-flipping.")
                        # DO NOT flip BETA - gwaslab already did it
                    } else {
                        # ALLELE_FLIPPED=FALSE or all NA, safe to flip downstream
                        if ("BETA.qtl" %in% names(merged_flip)) {
                            merged_flip$BETA.qtl <- -merged_flip$BETA.qtl
                        }
                        if ("BETA.gwas" %in% names(merged_flip)) {
                            merged_flip$BETA.gwas <- -merged_flip$BETA.gwas
                        }
                        merged_flip$downstream_flipped <- TRUE
                    }
                } else {
                    # No gwaslab tracking column - safe to apply downstream flip
                    if ("BETA.qtl" %in% names(merged_flip)) {
                        merged_flip$BETA.qtl <- -merged_flip$BETA.qtl
                    }
                    if ("BETA.gwas" %in% names(merged_flip)) {
                        merged_flip$BETA.gwas <- -merged_flip$BETA.gwas
                    }
                    merged_flip$downstream_flipped <- TRUE
                }
            }
            if (nrow(merged_forward) > 0) {
                merged_forward$allele_flipped <- FALSE
            }

            merged_dt <- rbindlist(list(merged_forward, merged_flip), fill = TRUE, use.names = TRUE)
             if (nrow(merged_dt) > 0) {
                 msg(glue("[MERGE] ✓ Strategy 3: {nrow(merged_dt)} SNPs matched (Forward: {nrow(merged_forward)}, Flipped: {nrow(merged_flip)})"))
                 merge_strategy <- "pos_allele"
                 # Priority: use GWAS rsID if available, otherwise use position-based ID
                 # This is critical for SuSiE which needs rsIDs to query PLINK
                 if ("rsid.gwas" %in% names(merged_dt)) {
                     merged_dt$snp <- merged_dt$rsid.gwas
                 } else if ("SNPID.gwas" %in% names(merged_dt)) {
                     merged_dt$snp <- merged_dt$SNPID.gwas
                 } else if ("SNP.gwas" %in% names(merged_dt)) {
                     merged_dt$snp <- merged_dt$SNP.gwas
                 } else if ("SNPID" %in% names(merged_dt)) {
                     merged_dt$snp <- merged_dt$SNPID
                 } else {
                     merged_dt$snp <- merged_dt$variant_id
                 }
                 snp_sample <- head(merged_dt$snp[!is.na(merged_dt$snp)], 3)
                 if (length(snp_sample) > 0 && !all(grepl("^rs", snp_sample))) {
                     msg(glue("[MERGE]   Warning: Non-rsID format detected: {paste(snp_sample, collapse=', ')}"))
                 }
             } else {
                 msg("[MERGE] ✗ Strategy 3 failed: No position + allele matches")
             }
        }

    }
    if (is.null(merged_dt) || nrow(merged_dt) == 0) {
        stop("*** MERGE FAILED: No overlapping SNPs found using any strategy! ***")
    }

    # =============================================================================
    # Final Summary
    # =============================================================================
    msg(glue("[MERGE] Final dataset: {nrow(merged_dt)} SNPs ready for colocalization"))

    p_gwas <- paste0(gwas_cols$pval, ".gwas")
    p_qtl  <- paste0(qtl_cols$pval, ".qtl")

    for (p_col in c(p_gwas, p_qtl)) {
        if (p_col %in% names(merged_dt)) {
            p_values <- merged_dt[[p_col]]
            needs_floor <- is.na(p_values) | p_values <= 0 | p_values < pvalue_floor
            n_fix <- sum(needs_floor, na.rm = TRUE)
            if (n_fix > 0) {
                merged_dt[needs_floor, (p_col) := pvalue_floor]
                msg(glue("  Replaced {n_fix} P values (NA/<=0/<floor) with {pvalue_floor} in {p_col}"))
            }
        }
    }
 
    # Fix column names that might have double suffixes
    # Find and rename columns like "P.gwas.gwas" to "P.gwas"
    double_suffix_cols <- grep("\\.gwas\\.gwas$|\\.qtl\\.qtl$", names(merged_dt), value = TRUE)
    if (length(double_suffix_cols) > 0) {
        msg(glue("  Fixing {length(double_suffix_cols)} double-suffix columns..."))
        for (col in double_suffix_cols) {
            new_col <- gsub("\\.gwas\\.gwas$", ".gwas", col)
            new_col <- gsub("\\.qtl\\.qtl$", ".qtl", new_col)
            setnames(merged_dt, col, new_col)
        }
    }
 
    # Find actual P column (GWAS P should be "P", QTL P should be "P.qtl")
    p_gwas_col <- if ("P" %in% names(merged_dt)) "P" else grep("^P\\.gwas$", names(merged_dt), value = TRUE)[1]
    if (is.na(p_gwas_col)) {
        p_gwas_col <- "P"  # Default fallback
    }
    msg(glue("  Using {p_gwas_col} for deduplication"))

    # Deduplication using data.table
    if (anyDuplicated(merged_dt$snp)) {
        n_dup <- sum(duplicated(merged_dt$snp))
        msg(glue("  Removing {n_dup} duplicate SNPs (keeping lowest {p_gwas_col})..."))

        # Use data.table for deduplication
        setorderv(merged_dt, p_gwas_col)
        merged_dt <- merged_dt[!duplicated(merged_dt$snp)]
    }

    # Quality Check
    required_cols <- c("snp",
                       paste0(c("BETA", "SE", "P"), ".gwas"),
                       paste0(c("BETA", "SE", "P"), ".qtl"))

    missing_cols <- setdiff(required_cols, names(merged_dt))
    if (length(missing_cols) > 0) {
        warning(glue("Missing columns: {paste(missing_cols, collapse=', ')}"))
    }

    msg(glue("Final dataset ready: {nrow(merged_dt)} SNPs\n"))
    return(merged_dt)
}

query_tabix_region <- function(tabix_file, chrom, start, end, col_names = NULL) {
  if (!file.exists(tabix_file)) {
      warning(glue("Tabix file not found: {tabix_file}"))
      return(NULL)
  }
  chrom_clean <- gsub("^chr", "", as.character(chrom))

  run_query <- function(query_chrom) {
      cmd <- glue("tabix {tabix_file} {query_chrom}:{start}-{end}")
      dt <- tryCatch({ suppressWarnings(fread(cmd = cmd, header = FALSE, showProgress = FALSE)) }, error = function(e) return(NULL))
      if (!is.null(dt) && nrow(dt) > 0) return(dt)
      return(NULL)
  }
  dt <- run_query(paste0("chr", chrom_clean))
  if (is.null(dt)) dt <- run_query(chrom_clean)
  return(dt)
}

get_ld_matrix <- function(variants, bfile, plink_bin, keep_file = NULL) {
    if (length(variants) == 0) return(NULL)
    fn <- tempfile(); on.exit(unlink(paste0(fn, "*")), add = TRUE)
    tryCatch({
        fwrite(list(variants), file = fn, col.names = FALSE)
        keep_cmd <- if (!is.null(keep_file) && file.exists(keep_file)) paste0("--keep ", keep_file) else ""
        cmd <- paste(plink_bin, "--bfile", bfile, "--extract", fn, keep_cmd, "--r square --make-just-bim --out", fn, "--silent")
        system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
        if (!file.exists(paste0(fn, ".ld"))) return(NULL)
        ld_mat <- as.matrix(fread(paste0(fn, ".ld")))
        bim <- fread(paste0(fn, ".bim"), col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"))
        colnames(ld_mat) <- bim$SNP; rownames(ld_mat) <- bim$SNP
        return(list(R = ld_mat, snp_info = bim))
    }, error = function(e) return(NULL))
}

write_coloc_result <- function(res_dt, output_file) {
  fwrite(res_dt, output_file, sep = "\t", quote = FALSE)
}

# ------------------------------------------------------------------------------
# merge_all_results: Aggregate colocalization results across all loci
# ------------------------------------------------------------------------------
merge_all_results <- function(output_dir,
                              pp4_threshold = 0.8,
                              merge_susie = TRUE,
                              save_summary = TRUE) {
    
    abf_dir <- file.path(output_dir, "abf")
    susie_dir <- file.path(output_dir, "susie")
    
    # Find all locus result files
    abf_files <- list.files(abf_dir, pattern = "_locus_results\\.csv$", full.names = TRUE)
    
    if (length(abf_files) == 0) {
        message("[SUM] No ABF results found to merge.")
        return(invisible(NULL))
    }
    
    message(glue("[SUM] Found {length(abf_files)} locus result file(s) to merge"))
    
    # Read and combine all ABF results
    all_results <- rbindlist(lapply(abf_files, function(f) {
        tryCatch({
            fread(f)
        }, error = function(e) {
            warning(glue("Failed to read {basename(f)}: {e$message}"))
            return(NULL)
        })
    }), fill = TRUE)
    
    if (nrow(all_results) == 0) {
        message("[SUM] No valid results found after merging.")
        return(invisible(NULL))
    }
    
    # Sort by PP4 (descending)
    all_results <- all_results[order(-PP4)]
    
    # Filter by PP4 threshold
    sig_results <- all_results[PP4 >= pp4_threshold]
    
    # Save merged results
    merged_file <- file.path(output_dir, "all_colocalization_results.csv")
    fwrite(all_results, merged_file)
    message(glue("[SUM] Saved merged results: {basename(merged_file)} ({nrow(all_results)} rows)"))
    
    # Save significant results
    if (nrow(sig_results) > 0) {
        sig_file <- file.path(output_dir, glue("significant_colocalizations_PP4_{pp4_threshold}.csv"))
        fwrite(sig_results, sig_file)
        message(glue("[SUM] Saved significant results (PP4 >= {pp4_threshold}): {basename(sig_file)} ({nrow(sig_results)} rows)"))
    } else {
        message(glue("[SUM] No significant results found with PP4 >= {pp4_threshold}"))
    }
    
    # Merge SuSiE results if requested
    if (merge_susie && dir.exists(susie_dir)) {
        susie_files <- list.files(susie_dir, pattern = "_susie\\.csv$", full.names = TRUE)
        
        if (length(susie_files) > 0) {
            message(glue("[SUM] Found {length(susie_files)} SuSiE result file(s)"))
            
            susie_results <- rbindlist(lapply(susie_files, function(f) {
                tryCatch({
                    dt <- fread(f)
                    dt$source_file <- basename(f)
                    return(dt)
                }, error = function(e) {
                    warning(glue("Failed to read {basename(f)}: {e$message}"))
                    return(NULL)
                })
            }), fill = TRUE)
            
            if (nrow(susie_results) > 0) {
                susie_merged_file <- file.path(output_dir, "all_susie_results.csv")
                fwrite(susie_results, susie_merged_file)
                message(glue("[SUM] Saved merged SuSiE results: {basename(susie_merged_file)} ({nrow(susie_results)} rows)"))
            }
        }
    }
    
    # Generate summary statistics if requested
    if (save_summary) {
        summary_stats <- list(
            total_tests = nrow(all_results),
            significant_colocalizations = nrow(sig_results),
            unique_gwas_traits = length(unique(all_results$GWAS_ID)),
            unique_qtl_datasets = length(unique(all_results$QTL_ID)),
            unique_loci = length(unique(all_results$Locus)),
            unique_genes = length(unique(all_results$Gene)),
            mean_pp4 = mean(all_results$PP4, na.rm = TRUE),
            median_pp4 = median(all_results$PP4, na.rm = TRUE),
            max_pp4 = max(all_results$PP4, na.rm = TRUE),
            mean_n_snps = mean(all_results$n_snps, na.rm = TRUE)
        )
        
        # Convert to data.frame for pretty printing
        summary_df <- data.frame(
            Metric = names(summary_stats),
            Value = unlist(summary_stats)
        )
        
        summary_file <- file.path(output_dir, "summary_statistics.txt")
        
        # Write formatted summary
        sink(summary_file)
        cat("=============================================================\n")
        cat("EasyColoc Analysis Summary\n")
        cat("=============================================================\n")
        cat(glue("Analysis completed: {Sys.time()}\n\n"))
        cat("Summary Statistics:\n")
        cat("-------------------------------------------------------------\n")
        print(summary_df, row.names = FALSE)
        cat("\n")
        
        # Top genes by PP4
        if (nrow(sig_results) > 0) {
            cat("\nTop 10 Gene-Locus Pairs by PP4:\n")
            cat("-------------------------------------------------------------\n")
            top_genes <- head(sig_results[, .(GWAS_ID, QTL_ID, Locus, Gene, PP4, n_snps)], 10)
            print(top_genes, row.names = FALSE)
        }
        
        cat("\n=============================================================\n")
        sink()
        
        message(glue("[SUM] Saved summary statistics: {basename(summary_file)}"))
    }
    
    message("[SUM] Merge complete!")
    return(invisible(all_results))
}
