# ------------------------------------------------------------------------------
# src/utils_analysis.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(dplyr)
  library(susieR)
})

find_lead_snp <- function(df, p_col = "P.gwas", snp_col = "snp") {
    if (!p_col %in% names(df)) stop(paste("Column", p_col, "not found"))
    df %>% arrange(.data[[p_col]]) %>% slice(1)
}

find_best_lead_snp_in_ld <- function(df, p_col = "P.gwas", snp_col = "snp", plink_bfile, plink_bin = "plink") {
    lead_row <- find_lead_snp(df, p_col, snp_col); lead_snp <- lead_row[[snp_col]]
    if (!exists("get_ld_matrix")) return(lead_snp)
    ld_check <- get_ld_matrix(lead_snp, plink_bfile, plink_bin)
    if (!is.null(ld_check) && nrow(ld_check$R) > 0) return(lead_snp) 
    top_candidates <- df %>% arrange(.data[[p_col]]) %>% slice(1:100) %>% pull(.data[[snp_col]])
    ld_proxy <- get_ld_matrix(top_candidates, plink_bfile, plink_bin)
    if (is.null(ld_proxy) || nrow(ld_proxy$R) == 0) return(lead_snp) 
    available <- colnames(ld_proxy$R)
    best <- df %>% filter(.data[[snp_col]] %in% available) %>% arrange(.data[[p_col]]) %>% slice(1) %>% pull(.data[[snp_col]])
    return(best)
}

get_coloc_results <- function(colocInputFile, gwas_type = "cc", gwas_prop = 0.5, gwas_N = NULL, qtl_N = NULL, use_susie = TRUE, susie_threshold = 0.75, plink_bfile = NULL, plink_bin = "plink") { 
    get_clean_maf <- function(df, suffix) {
        candidates <- c(paste0("EAF", suffix), paste0("MAF", suffix), paste0("af", suffix))
        if (suffix == ".gwas") candidates <- c(candidates, "EAF", "MAF", "af")
        raw_freq <- NULL
        for (col in candidates) if (col %in% names(df)) { raw_freq <- df[[col]]; break }
        if (is.null(raw_freq)) return(rep(0.1, nrow(df)))
        freq <- as.numeric(raw_freq); maf <- pmin(freq, 1 - freq)
        if (any(is.na(maf))) maf[is.na(maf)] <- 0.05
        epsilon <- 1e-6; maf[maf < epsilon] <- epsilon; maf[maf > (1 - epsilon)] <- 1 - epsilon
        return(maf)
    }

    N_gwas_val <- if("N.gwas" %in% names(colocInputFile)) colocInputFile$N.gwas else if ("N" %in% names(colocInputFile)) colocInputFile$N else gwas_N
    if (is.null(N_gwas_val)) stop("GWAS N missing.")
    N_qtl_val <- if("N.qtl" %in% names(colocInputFile)) colocInputFile$N.qtl else qtl_N
    if (is.null(N_qtl_val)) stop("QTL N missing.")
    
    d1 <- list(snp = as.character(colocInputFile$snp), beta = as.numeric(colocInputFile$BETA.gwas), varbeta = as.numeric(colocInputFile$SE.gwas)^2, pvalues = as.numeric(colocInputFile$P.gwas), type = gwas_type, N = as.numeric(N_gwas_val), MAF = get_clean_maf(colocInputFile, ".gwas"))
    if (gwas_type == "cc") { if (is.null(gwas_prop)) stop("Case prop required."); d1$s <- as.numeric(gwas_prop) }
    d2 <- list(snp = as.character(colocInputFile$snp), beta = as.numeric(colocInputFile$BETA.qtl), varbeta = as.numeric(colocInputFile$SE.qtl)^2, pvalues = as.numeric(colocInputFile$P.qtl), type = "quant", N = as.numeric(N_qtl_val), MAF = get_clean_maf(colocInputFile, ".qtl"))
    
    message("Running Coloc ABF...")
    res_abf <- coloc.abf(d1, d2)
    pp4 <- as.numeric(res_abf$summary["PP.H4.abf"])
    
    if (!use_susie || is.na(pp4) || pp4 < susie_threshold) {
        return(res_abf)
    }
    
    message(sprintf("  -> Strong signal detected (PP.H4 = %.4f). Running SuSiE...", pp4))
    
    if (!exists("get_ld_matrix") || is.null(plink_bfile)) { warning("LD tools missing. Skipping SuSiE."); return(res_abf) }
     
    ld_res <- get_ld_matrix(d1$snp, plink_bfile, plink_bin)
    
    if (is.null(ld_res) || is.null(ld_res$R)) { 
        warning(glue("LD calculation failed (NULL matrix). Check if SNP IDs (e.g. {d1$snp[1]}) match PLINK reference."))
        return(res_abf) 
    }
    
    common_snps <- intersect(d1$snp, colnames(ld_res$R))
    if (length(common_snps) < 5) { 
        warning(glue("Too few overlapping SNPs ({length(common_snps)}) between data and LD ref. Skipping SuSiE."))
        return(res_abf) 
    }
    
    filter_data <- function(d, snps) {
        idx <- match(snps, d$snp)
        d$snp <- d$snp[idx]; d$beta <- d$beta[idx]; d$varbeta <- d$varbeta[idx]
        if(!is.null(d$MAF)) d$MAF <- d$MAF[idx]
        return(d)
    }
    
    d1_susie <- filter_data(d1, common_snps); d2_susie <- filter_data(d2, common_snps)
    R_susie <- ld_res$R[common_snps, common_snps]
    
    susie_1 <- try(susie_rss(bhat = d1_susie$beta, shat = sqrt(d1_susie$varbeta), R = R_susie, n = d1_susie$N[1], verbose = FALSE), silent = TRUE)
    susie_2 <- try(susie_rss(bhat = d2_susie$beta, shat = sqrt(d2_susie$varbeta), R = R_susie, n = d2_susie$N[1], verbose = FALSE), silent = TRUE)
    
    if (inherits(susie_1, "try-error") || inherits(susie_2, "try-error")) { warning("SuSiE failed to converge."); return(res_abf) }
    
    res_susie <- coloc.susie(susie_1, susie_2)
    message("  SuSiE Coloc complete.")
    
    res_abf$susie_result <- res_susie
    return(res_abf)
}
