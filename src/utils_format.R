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
                                      dataset_id = NULL) {

  message("--- Starting Harmonization & LiftOver ---")

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
    print("Loading sumstats with GWASLab...")
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
        print("Running GWASLab harmonization...")
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
        print("Performing LiftOver from " + src_b + " to " + tgt_b + "...")
        mysumstats.liftover(from_build=src_b, to_build=tgt_b, remove=True)

    mysumstats.data.to_csv("%s", sep="\\t", index=False, float_format="%%.6g")

except Exception as e:
    print(f"Error in GWASLab: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
', input_file, source_build,
   ref_fasta, dbsnp_pattern, vcf_pattern, ref_alt_freq, source_build, target_build,
   n_threads,
   output_file)

     writeLines(py_script, py_script_file)

  cmd <- glue("micromamba run -n {env_name} python {py_script_file}")
  message("Running Python script...")
  exit_code <- system(cmd)

   if (exit_code == 0) {
       if (file.exists(output_file)) {
         res_dt <- fread(output_file)
         message(glue("Processing complete. Input SNPs: {nrow(sumstats_dt)} -> Output SNPs: {nrow(res_dt)}"))
         message(glue("  Output columns: {paste(names(res_dt), collapse=', ')}"))
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
                                  use_hash_table = TRUE,
                                  min_snps = 30,
                                  pvalue_floor = 1e-300) {

    message("\n--- MERGING DATA ---")
    gwas_df <- as.data.table(gwas_df)
    qtl_df <- as.data.table(qtl_df)
    message(glue("Input: GWAS={nrow(gwas_df)} SNPs, QTL={nrow(qtl_df)} SNPs"))

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
         qtl_df[, `:=`(EA = ALT, NEA = REF)]
    }
    merged_dt <- NULL
    merge_strategy <- "none"

    if (gwas_has_rsid && qtl_has_rsid) {
        message("Match SNPs...")

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
        message("  Matching SNPs...")
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

            message(glue("rsID match: {nrow(merged_dt)} SNPs"))
            merge_strategy <- "rsid_direct"
            if ("SNPID.gwas" %in% names(merged_dt)) {
                merged_dt$snp <- merged_dt$SNPID.gwas
            } else if ("SNP.gwas" %in% names(merged_dt)) {
                merged_dt$snp <- merged_dt$SNP.gwas
            } else if ("SNPID" %in% names(merged_dt)) {
                merged_dt$snp <- merged_dt$SNPID
            } else if ("SNP" %in% names(merged_dt)) {
                merged_dt$snp <- merged_dt$SNP
            } else {
                merged_dt$snp <- merged_dt[[gwas_id_col]]
            }
        } else {
            message("No rsID overlap detected")
            merged_dt <- NULL
        }
    }

    if ((is.null(merged_dt) || nrow(merged_dt) < min_snps) && use_hash_table && exists("convert_rsid_to_pos")) {
        message("Attempting rsID -> chr_pos conversion...")
        gwas_id_col <- if ("rsid" %in% names(gwas_df)) "rsid" else if ("SNPID" %in% names(gwas_df)) "SNPID" else if ("SNP" %in% names(gwas_df)) "SNP" else NULL

        if (!is.null(gwas_id_col) && gwas_has_pos) {
            message("  Converting GWAS rsIDs to chr_pos format...")
            
            gwas_chrpos <- convert_rsid_to_pos(
                gwas_df[[gwas_id_col]],
                gwas_df$CHR
            )
            gwas_df$gwas_chrpos_id <- gwas_chrpos

            n_converted <- sum(!is.na(gwas_chrpos))
            message(glue("  Converted {n_converted}/{nrow(gwas_df)} GWAS SNPs ({round(n_converted/nrow(gwas_df)*100, 1)}%)"))
            
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
                
                 if (nrow(merged_dt) > 0) {
                     message(glue("  ✓ Hash table conversion match: {nrow(merged_dt)} SNPs"))
                     merge_strategy <- "hash_conversion"
                     rsid_col <- if ("SNPID.gwas" %in% names(merged_dt)) "SNPID.gwas"
                                 else if ("SNP.gwas" %in% names(merged_dt)) "SNP.gwas"
                                 else if ("SNPID" %in% names(merged_dt)) "SNPID"
                                 else if ("SNP" %in% names(merged_dt)) "SNP"
                                 else NULL
                     
                     if (!is.null(rsid_col)) {
                         merged_dt$snp <- merged_dt[[rsid_col]]
                         message(glue("  Using {rsid_col} for SNP column"))
                     } else {
                         setnames(merged_dt, "gwas_chrpos_id", "snp")
                     }
                     
                     snp_sample <- head(merged_dt$snp[!is.na(merged_dt$snp)], 3)
                     if (length(snp_sample) > 0 && !all(grepl("^rs", snp_sample))) {
                         message(glue("  [WARN] SNP not in rsID format: {paste(snp_sample, collapse=', ')}"))
                     }
                 }
            }
        }
    }
    if ((is.null(merged_dt) || nrow(merged_dt) < min_snps) && gwas_has_pos && qtl_has_pos) {
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
            merged_forward <- inner_join(
                gwas_df, qtl_df,
                by = "variant_id",
                suffix = c(".gwas", ".qtl"),
                relationship = "many-to-many"
            )
            merged_forward <- as.data.table(merged_forward)

            # Flip Match
            merged_flip <- inner_join(
                gwas_df, qtl_df,
                by = c("variant_id_flip" = "variant_id"),
                suffix = c(".gwas", ".qtl"),
                relationship = "many-to-many"
            )
            merged_flip <- as.data.table(merged_flip)

            # Mark flipped SNPs
            if (nrow(merged_flip) > 0) {
                merged_flip$allele_flipped <- TRUE
                if ("BETA.qtl" %in% names(merged_flip)) {
                    merged_flip$BETA.qtl <- -merged_flip$BETA.qtl
                }
                if ("BETA.gwas" %in% names(merged_flip)) {
                    merged_flip$BETA.gwas <- -merged_flip$BETA.gwas
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
                 if ("SNPID.gwas" %in% names(merged_dt)) {
                     merged_dt$snp <- merged_dt$SNPID.gwas
                     message("  Using SNPID.gwas for SNP column")
                 } else if ("SNP.gwas" %in% names(merged_dt)) {
                     merged_dt$snp <- merged_dt$SNP.gwas
                     message("  Using SNP.gwas for SNP column")
                 } else if ("rsid.gwas" %in% names(merged_dt)) {
                     merged_dt$snp <- merged_dt$rsid.gwas
                     message("  Using rsid.gwas for SNP column")
                 } else {
                     merged_dt$snp <- merged_dt$variant_id
                     message("  No rsID column found, using variant_id")
                 }
                 snp_sample <- head(merged_dt$snp[!is.na(merged_dt$snp)], 3)
                 if (length(snp_sample) > 0 && !all(grepl("^rs", snp_sample))) {
                     message(glue("  [WARN] SNP not in rsID format: {paste(snp_sample, collapse=', ')}"))
                 }
             }
        }
    }

    if (is.null(merged_dt) || nrow(merged_dt) == 0) {
        stop("*** MERGE FAILED: No overlapping SNPs found using any strategy! ***")
    }

    message(glue("Merged dataset: {nrow(merged_dt)} SNPs"))
    message(glue("  SNP column sample: {paste(head(merged_dt$snp, 5), collapse=', ')}"))

    p_gwas <- paste0(gwas_cols$pval, ".gwas")
    p_qtl  <- paste0(qtl_cols$pval, ".qtl")

    for (p_col in c(p_gwas, p_qtl)) {
        if (p_col %in% names(merged_dt)) {
            p_values <- merged_dt[[p_col]]
            needs_floor <- is.na(p_values) | p_values <= 0 | p_values < pvalue_floor
            n_fix <- sum(needs_floor, na.rm = TRUE)
            if (n_fix > 0) {
                merged_dt[needs_floor, (p_col) := pvalue_floor]
                message(glue("  Replaced {n_fix} P values (NA/<=0/<floor) with {pvalue_floor} in {p_col}"))
            }
        }
    }
 
    # Fix column names that might have double suffixes
    # Find and rename columns like "P.gwas.gwas" to "P.gwas"
    double_suffix_cols <- grep("\\.gwas\\.gwas$|\\.qtl\\.qtl$", names(merged_dt), value = TRUE)
    if (length(double_suffix_cols) > 0) {
        message(glue("  Fixing {length(double_suffix_cols)} double-suffix columns..."))
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
    message(glue("  Using {p_gwas_col} for deduplication"))

    # Deduplication using data.table
    if (anyDuplicated(merged_dt$snp)) {
        n_dup <- sum(duplicated(merged_dt$snp))
        message(glue("  Removing {n_dup} duplicate SNPs (keeping lowest {p_gwas_col})..."))

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
