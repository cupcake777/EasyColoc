# ------------------------------------------------------------------------------
# src/utils_helpers.R
# Des: print info/liftover/recom
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(glue)
  library(data.table)
  library(dplyr)
})

print_config_settings <- function() {
    if (getOption("coloc_config_printed", default = FALSE)) {
        return(invisible(NULL))
    }
    get_var <- function(var_name) {
        if (exists(var_name, envir = .GlobalEnv)) {
            val <- get(var_name, envir = .GlobalEnv)
            if (is.null(val)) return("NULL")
            return(as.character(val))
        }
        return("N/A")
    }

    cat(glue("
    ========================================================
    [CONFIG SUMMARY] {Sys.time()}
    ========================================================
    Trait Info:
      - Name: {get_var('trait')}
      - File: {basename(get_var('traitFilePath'))}
      - Type: {get_var('traitType')} (Prop: {get_var('traitProp')})

    Column Mapping (Standardized):
      - SNP:  {get_var('trait_SNPcol')}
      - Pval: {get_var('trait_Pcol')}
      - Beta: {get_var('trait_BETAcol')}
      - N:    {get_var('trait_Ncol')}

    Locus Info:
      - Region: Chr {get_var('chrom')}:{get_var('colocStart')} - {get_var('colocStop')}
      - Lead SNP: {get_var('lead_SNP')}
      - Build: {get_var('trait_build')}

    System:
      - LiftOver Chain: {basename(get_var('liftOver_chain'))}
    ========================================================
    \n"))
    options(coloc_config_printed = TRUE)
}

get_liftover_point <- function(chrom, pos, chain_file) {
    if (is.null(chain_file) || !file.exists(chain_file)) stop("LiftOver chain file not found.")

    chrom_str <- if (!grepl("^chr", chrom)) paste0("chr", chrom) else chrom
    
    bed_df <- data.frame(
        chr = chrom_str,
        start = as.integer(pos) - 1,
        end = as.integer(pos),
        id = "pt"
    )
    
    f_in <- tempfile(fileext = ".bed")
    f_out <- tempfile(fileext = ".bed")
    f_unmapped <- tempfile(fileext = ".unmapped")
    
    on.exit({
        if (file.exists(f_in)) unlink(f_in)
        if (file.exists(f_out)) unlink(f_out)
        if (file.exists(f_unmapped)) unlink(f_unmapped)
    })
    
    fwrite(bed_df, f_in, sep = "\t", col.names = FALSE, scipen = 50)
    cmd <- glue("liftOver {f_in} {chain_file} {f_out} {f_unmapped}")
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    if (file.size(f_out) == 0) return(NULL)
    
    res_df <- fread(f_out, header = FALSE, col.names = c("chr", "start", "end", "id"))
    
    return(list(
        chrom = gsub("^chr", "", res_df$chr[1]),
        pos = res_df$end[1]
    ))
}

load_recomb_map <- function(chrom, start, end, recomb_prefix) {
    if (is.null(recomb_prefix)) return(NULL)
    chrom_clean <- gsub("chr", "", as.character(chrom))
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
                return(dt %>% filter(`Position(bp)` >= start & `Position(bp)` <= end))
            }
        }
    }
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
                    select(`Position(bp)`, `Rate(cM/Mb)`) %>%
                    filter(`Position(bp)` >= start & `Position(bp)` <= end)
                
                return(dt_out)
            }
        }
    }
    return(NULL)
}
