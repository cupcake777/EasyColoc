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
    maf  = "EAF"
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

  # Now safe to use := because df_std is guaranteed to be a data.table
  if (case_control && "BETA" %notin% names(df_std) && "OR" %in% names(df_std)) {
    df_std[, BETA := log(OR)]
  }

  return(df_std)
}

run_gwaslab_harmonization <- function(sumstats_dt,
                                      ref_fasta,
                                      source_build = "19",
                                      target_build = "38",
                                      env_name = "gwaslab",
                                      save_dir = NULL,
                                      dataset_id = NULL) {

  message("--- Starting GWASLab Harmonization & LiftOver ---")

  if (!file.exists(ref_fasta)) {
    warning(glue("Reference FASTA not found: {ref_fasta}. Skipping harmonization."))
    return(sumstats_dt)
  }

  use_cache <- !is.null(save_dir) && !is.null(dataset_id)
  final_output_file <- NULL

  if (use_cache) {
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
      ref_name <- tools::file_path_sans_ext(basename(ref_fasta))
      # Cache filename includes build info
      final_output_file <- file.path(save_dir, glue("{dataset_id}_{ref_name}_b{source_build}to{target_build}_harmonized.tsv"))

      if (file.exists(final_output_file)) {
          message(glue("Found cached harmonized file: {final_output_file}"))
          message("Loading from cache...")
          return(fread(final_output_file))
      }
  }

  temp_id <- format(Sys.time(), "%H%M%S")
  input_file <- glue("temp_gwas_in_{temp_id}.tsv")
  output_file <- if(use_cache) final_output_file else glue("temp_gwas_out_{temp_id}.tsv")
  py_script_file <- glue("run_gwaslab_{temp_id}.py")

  on.exit({
    if (file.exists(input_file)) file.remove(input_file)
    if (file.exists(py_script_file)) file.remove(py_script_file)
    if (!use_cache && file.exists(output_file)) file.remove(output_file)
  }, add = TRUE)

  fwrite(sumstats_dt, input_file, sep = "\t", na = "NA", quote = FALSE)

  # Python script with conditional LiftOver
  py_code <- glue("
import gwaslab as gl
import pandas as pd
import sys
try:
    print('Loading sumstats with GWASLab...')
    # Load assuming gwaslab format since we formatted it previously
    # Explicitly set the source build
    mysumstats = gl.Sumstats('{input_file}', fmt='gwaslab', build='{source_build}')
    mysumstats.basic_check()
    
    print(f'Harmonizing against: {ref_fasta}')
    mysumstats.harmonize(ref_seq='{ref_fasta}', ref_infer=True, remove=True)
    
    # Conditional LiftOver
    src_b = '{source_build}'
    tgt_b = '{target_build}'
    
    if src_b != tgt_b:
        print(f'Performing LiftOver from {{src_b}} to {{tgt_b}}...')
        mysumstats.liftover(n_cores=3, from_build=src_b, to_build=tgt_b, remove=True)
    
    mysumstats.data.to_csv('{output_file}', sep='\\t', index=False, float_format='%.6g')

except Exception as e:
    print(f'Error in GWASLab: {{e}}')
    sys.exit(1)
")

  writeLines(py_code, py_script_file)

  cmd <- glue("micromamba run -n {env_name} python {py_script_file}")
  message("Running Python script...")
  exit_code <- system(cmd)

  if (exit_code == 0) {
      if (file.exists(output_file)) {
        res_dt <- fread(output_file)
        message(glue("Processing complete. Input SNPs: {nrow(sumstats_dt)} -> Output SNPs: {nrow(res_dt)}"))
        if(use_cache) message(glue("Saved result to: {output_file}"))
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
                                  use_hash_table = TRUE) {

    message("\n=== MERGING GWAS AND QTL DATA ===")
    gwas_df <- as.data.table(gwas_df)
    qtl_df <- as.data.table(qtl_df)
    message(glue("Input: GWAS={nrow(gwas_df)} SNPs, QTL={nrow(qtl_df)} SNPs"))
    gwas_has_rsid <- any(c("SNPID", "SNP", "rsid") %in% names(gwas_df))
    qtl_has_rsid <- any(c("SNPID", "variant_id", "snp") %in% names(qtl_df))

    gwas_has_pos <- all(c("CHR", "POS") %in% names(gwas_df))
    qtl_has_pos <- all(c("CHR", "POS") %in% names(qtl_df)) ||
                   all(c("chr", "end") %in% names(qtl_df))
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
         qtl_df[, `:=`(EA = ALT, NEA = REF)] # Assuming ALT is Effect Allele for QTL usually
    }
    merged_dt <- NULL
    merge_strategy <- "none"

    if (gwas_has_rsid && qtl_has_rsid) {
        message("Match SNPs...")

        gwas_id_col <- if ("SNPID" %in% names(gwas_df)) "SNPID" else if ("SNP" %in% names(gwas_df)) "SNP" else "rsid"
        qtl_id_col <- if ("SNPID" %in% names(qtl_df)) "SNPID" else if ("variant_id" %in% names(qtl_df)) "variant_id" else "snp"

        gwas_sample <- head(gwas_df[[gwas_id_col]], 3)
        qtl_sample <- head(qtl_df[[qtl_id_col]], 3)
        message(glue("  GWAS IDs: {paste(gwas_sample, collapse=', ')}"))
        message(glue("  QTL IDs:  {paste(qtl_sample, collapse=', ')}"))

        merged_dt <- merge(
            gwas_df, qtl_df,
            by.x = gwas_id_col,
            by.y = qtl_id_col,
            suffixes = c(".gwas", ".qtl")
        )

        if (nrow(merged_dt) > 0) {
            message(glue("Direct rsID match: {nrow(merged_dt)} SNPs"))
            merge_strategy <- "rsid_direct"
            if (gwas_id_col %in% names(merged_dt)) setnames(merged_dt, gwas_id_col, "snp")
        } else {
            message("No direct rsID overlap detected")
            merged_dt <- NULL 
        }
    }

    if ((is.null(merged_dt) || nrow(merged_dt) < 30) && use_hash_table && exists("convert_rsid_to_chrpos_optimized")) {
        message("Attempting rsID -> chr_pos conversion...")
        gwas_id_col <- if ("SNPID" %in% names(gwas_df)) "SNPID" else if ("SNP" %in% names(gwas_df)) "SNP" else if ("rsid" %in% names(gwas_df)) "rsid" else NULL

        if (!is.null(gwas_id_col) && gwas_has_pos) {
            message("  Converting GWAS rsIDs to chr_pos format (chromosome-aware)...")
            
            gwas_chrpos <- convert_rsid_to_chrpos_optimized(
                gwas_df[[gwas_id_col]],
                gwas_df$CHR
            )
            gwas_df$gwas_chrpos_id <- gwas_chrpos # Works on DT and DF

            n_converted <- sum(!is.na(gwas_chrpos))
            message(glue("  Converted {n_converted}/{nrow(gwas_df)} GWAS SNPs ({round(n_converted/nrow(gwas_df)*100, 1)}%)"))
            
            qtl_id_col <- if ("variant_id" %in% names(qtl_df)) "variant_id" else if ("snp" %in% names(qtl_df)) "snp" else NULL

            if (!is.null(qtl_id_col)) {
                merged_dt <- merge(
                    gwas_df, qtl_df,
                    by.x = "gwas_chrpos_id",
                    by.y = qtl_id_col,
                    suffixes = c(".gwas", ".qtl")
                )
                
                if (nrow(merged_dt) > 0) {
                    message(glue("  ✓ Hash table conversion match: {nrow(merged_dt)} SNPs"))
                    merge_strategy <- "hash_conversion"
                    
                    if (gwas_id_col %in% names(gwas_df)) {
                         idx <- match(merged_dt$gwas_chrpos_id, gwas_df$gwas_chrpos_id)
                         merged_dt$snp <- gwas_df[[gwas_id_col]][idx]
                    } else {
                         setnames(merged_dt, "gwas_chrpos_id", "snp")
                    }
                }
            }
        }
    }
    if ((is.null(merged_dt) || nrow(merged_dt) < 30) && gwas_has_pos && qtl_has_pos) {
        message("\n[Strategy 3] Attempting position + allele matching...")
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
            merged_forward <- merge(
                gwas_df, qtl_df,
                by = "variant_id",
                suffixes = c(".gwas", ".qtl")
            )
            
            # Flip Match
            merged_flip <- merge(
                gwas_df, qtl_df,
                by.x = "variant_id_flip",
                by.y = "variant_id",
                suffixes = c(".gwas", ".qtl")
            )

            # Mark flipped SNPs
            if (nrow(merged_flip) > 0) {
                merged_flip$allele_flipped <- TRUE
                # Safe flip: check if column exists before flipping
                if ("BETA.qtl" %in% names(merged_flip)) {
                    merged_flip$BETA.qtl <- -merged_flip$BETA.qtl
                }
            }
            if (nrow(merged_forward) > 0) {
                merged_forward$allele_flipped <- FALSE
            }

            merged_dt <- rbindlist(list(merged_forward, merged_flip), fill = TRUE, use.names = TRUE)

            if (nrow(merged_dt) > 0) {
                message(glue("  ✓ Position+Allele match: {nrow(merged_dt)} SNPs"))
                message(glue("    (Forward: {nrow(merged_forward)}, Flipped: {nrow(merged_flip)})"))
                merge_strategy <- "pos_allele"
                
                # Assign 'snp' column
                if ("SNPID" %in% names(gwas_df)) {
                    # Try to map back to rsID using variant_id
                    # This assumes variant_id is unique enough for lookup or just use first match
                    idx <- match(merged_dt$variant_id, gwas_df$variant_id)
                    merged_dt$snp <- gwas_df$SNPID[idx]
                    merged_dt$snp[is.na(merged_dt$snp)] <- merged_dt$variant_id[is.na(merged_dt$snp)]
                } else {
                    merged_dt$snp <- merged_dt$variant_id
                }
            }
        }
    }

    if (is.null(merged_dt) || nrow(merged_dt) == 0) {
        stop("*** MERGE FAILED: No overlapping SNPs found using any strategy! ***")
    }

  message(glue("Merged dataset: {nrow(merged_dt)} SNPs"))
    p_gwas <- paste0(gwas_cols$pval, ".gwas")
    p_qtl  <- paste0(qtl_cols$pval, ".qtl")

    for (p_col in c(p_gwas, p_qtl)) {
        if (p_col %in% names(merged_dt)) {
            # Use data.table assignment for speed/safety
            merged_dt[get(p_col) == 0, (p_col) := 1e-300]
            n_fix <- sum(merged_dt[[p_col]] == 1e-300, na.rm=TRUE)
            if (n_fix > 0) message(glue("  Replaced {n_fix} P=0 values with 1e-300 in {p_col}"))
        }
    }

    # Deduplication
    if (anyDuplicated(merged_dt$snp)) {
        n_dup <- sum(duplicated(merged_dt$snp))
        message(glue("  Removing {n_dup} duplicate SNPs (keeping lowest P.gwas)..."))

        # Use dplyr syntax but ensure object is compatible
        merged_dt <- merged_dt %>%
            group_by(snp) %>%
            slice_min(order_by = !!sym(p_gwas), n = 1, with_ties = FALSE) %>%
            ungroup() %>%
            as.data.table()
    }

    # Quality Check
    required_cols <- c("snp",
                       paste0(c("BETA", "SE", "P"), ".gwas"),
                       paste0(c("BETA", "SE", "P"), ".qtl"))

    missing_cols <- setdiff(required_cols, names(merged_dt))
    if (length(missing_cols) > 0) {
        warning(glue("Missing columns: {paste(missing_cols, collapse=', ')}"))
    }

    message(glue("Final dataset ready: {nrow(merged_dt)} SNPs\n"))
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

get_ld_matrix <- function(variants, bfile, plink_bin) {
    if (length(variants) == 0) return(NULL)
    fn <- tempfile(); on.exit(unlink(paste0(fn, "*")), add = TRUE)
    tryCatch({
        fwrite(list(variants), file = fn, col.names = FALSE)
        cmd <- paste(plink_bin, "--bfile", bfile, "--extract", fn, "--r square --make-just-bim --out", fn, "--silent")
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
