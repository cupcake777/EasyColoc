# ------------------------------------------------------------------------------
# src/utils_helpers.R
# Des: filename sanitization, recombination map loading
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(glue)
  library(data.table)
  library(dplyr)
})

# Sanitize filename by removing problematic characters
sanitize_filename <- function(filename) {
  # Replace problematic characters with underscores
  sanitized <- gsub("[^[:alnum:]_.-]", "_", filename)
  # Replace multiple underscores with single underscore
  sanitized <- gsub("_+", "_", sanitized)
  # Remove leading/trailing underscores
  sanitized <- gsub("^_|_$", "", sanitized)
  return(sanitized)
}

.recomb_cache <- new.env(parent = emptyenv())
.plink_bim_cache <- new.env(parent = emptyenv())

load_recomb_map <- function(chrom, start, end, recomb_prefix) {
    if (is.null(recomb_prefix)) return(NULL)
    chrom_clean <- gsub("chr", "", as.character(chrom))
    
    # Check cache
    dt_full <- NULL
    if (exists(chrom_clean, envir = .recomb_cache)) {
        dt_full <- get(chrom_clean, envir = .recomb_cache)
    }
    
    if (is.null(dt_full)) {
        fname_txt <- paste0(recomb_prefix, "_recombination_map_hapmap_format_hg38_chr_", chrom_clean, ".txt")
        if (file.exists(fname_txt)) {
            dt <- tryCatch({ fread(fname_txt, header = TRUE) }, error = function(e) NULL)
            if (!is.null(dt) && nrow(dt) > 0) {
                if (!"Position(bp)" %in% names(dt)) {
                    if (ncol(dt) >= 3) {
                        names(dt)[2] <- "Position(bp)"
                        names(dt)[3] <- "Rate(cM/Mb)"
                    }
                }
                if ("Position(bp)" %in% names(dt)) {
                     dt_full <- dt
                }
            }
        }
        
        if (is.null(dt_full)) {
            fname_bed <- paste0(recomb_prefix, "_recombination_map_hg38_chr_", chrom_clean, ".bed")
            if (file.exists(fname_bed)) {
                dt <- tryCatch({ fread(fname_bed, header = TRUE) }, error = function(e) NULL)
                if (!is.null(dt) && nrow(dt) > 0) {
                    if (ncol(dt) >= 4) {
                        dt_out <- dt %>%
                            select(
                                `Position(bp)` = 2, 
                                raw_rate = 4
                            ) %>%
                            mutate(
                                `Rate(cM/Mb)` = raw_rate * 1e8
                            ) %>%
                            select(`Position(bp)`, `Rate(cM/Mb)`)
                        dt_full <- dt_out
                    }
                }
            }
        }
        
        # Cache and cleanup
        if (!is.null(dt_full)) {
             # Clear other chromosomes to save memory
             rm(list = setdiff(ls(.recomb_cache), chrom_clean), envir = .recomb_cache)
             assign(chrom_clean, dt_full, envir = .recomb_cache)
        }
    }

    if (!is.null(dt_full)) {
         return(dt_full)
    }
    
    return(NULL)
}

load_plink_bim_index <- function(plink_bfile) {
    if (is.null(plink_bfile) || !nzchar(plink_bfile)) {
        return(NULL)
    }
    cache_key <- normalizePath(plink_bfile, mustWork = FALSE)
    if (exists(cache_key, envir = .plink_bim_cache, inherits = FALSE)) {
        return(get(cache_key, envir = .plink_bim_cache))
    }

    bim_file <- paste0(plink_bfile, ".bim")
    if (!file.exists(bim_file)) {
        warning(glue("PLINK BIM file not found: {bim_file}"))
        return(NULL)
    }

    ref_panel <- tryCatch(
        fread(
            bim_file,
            header = FALSE,
            select = c(1, 2, 4, 5, 6),
            col.names = c("CHR", "plink_snp_id", "POS", "V5", "V6")
        ),
        error = function(e) NULL
    )

    if (is.null(ref_panel) || nrow(ref_panel) == 0) {
        warning(glue("Failed to load PLINK BIM index: {bim_file}"))
        return(NULL)
    }

    ref_panel[, CHR := gsub("^chr", "", as.character(CHR), ignore.case = TRUE)]
    ref_panel[, POS := suppressWarnings(as.numeric(POS))]
    assign(cache_key, ref_panel, envir = .plink_bim_cache)
    ref_panel
}
