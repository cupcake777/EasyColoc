# ------------------------------------------------------------------------------
# src/utils_analysis.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(dplyr)
  library(susieR)
})

get_coloc_results <- function(colocInputFile, gwas_type = "cc", gwas_prop = 0.5, gwas_N = NULL, qtl_N = NULL, use_susie = TRUE, susie_threshold = 0.75, plink_bfile = NULL, plink_bin = "plink", p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, maf_default = 0.1, maf_na_replacement = 0.05, maf_epsilon = 1e-6, keep_file = NULL) {
    # ==========================================================================
    # get_clean_maf: Extract or impute minor allele frequencies
    # ==========================================================================
    # Priority: (1) EAF/MAF/af column if present
    #           (2) Default maf_default (configurable, default 0.1)
    # Missing values replaced with maf_na_replacement (default 0.05)
    # Values clamped to [maf_epsilon, 1-maf_epsilon] for numerical stability
    get_clean_maf <- function(df, suffix) {
        candidates <- c(paste0("EAF", suffix), paste0("MAF", suffix), paste0("af", suffix))
        if (suffix == ".gwas") candidates <- c(candidates, "EAF", "MAF", "af")
        raw_freq <- NULL
        for (col in candidates) if (col %in% names(df)) { raw_freq <- df[[col]]; break }

        # Case 1: No frequency column found - use default
        if (is.null(raw_freq)) {
            warning(sprintf("[MAF] No allele frequency column found; using default MAF=%.2f for all %d SNPs",
                           maf_default, nrow(df)))
            return(rep(maf_default, nrow(df)))
        }

        freq <- as.numeric(raw_freq); maf <- pmin(freq, 1 - freq)

        # Case 2: Missing values in frequency column
        n_na_maf <- sum(is.na(maf))
        if (n_na_maf > 0) {
            message(sprintf("[MAF] %d SNPs (%.1f%%) had missing MAF - replaced with %.3f",
                           n_na_maf, n_na_maf/nrow(df)*100, maf_na_replacement))
            maf[is.na(maf)] <- maf_na_replacement
        }

        # Clamp to valid range to prevent numerical issues
        maf[maf < maf_epsilon] <- maf_epsilon
        maf[maf > (1 - maf_epsilon)] <- 1 - maf_epsilon
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
    res_abf <- coloc.abf(d1, d2, p1 = p1, p2 = p2, p12 = p12)
    pp4 <- as.numeric(res_abf$summary["PP.H4.abf"])
    
    if (!use_susie || is.na(pp4) || pp4 < susie_threshold) {
        return(res_abf)
    }
    
    message(sprintf("  -> Strong signal detected (PP.H4 = %.4f). Running SuSiE...", pp4))
    
    # Check if SNPs are in rsID format (required for PLINK)
    snp_sample <- head(d1$snp, 5)
    is_rsid_format <- all(grepl("^rs", snp_sample), na.rm = TRUE)
    if (!is_rsid_format) {
        message("SNPs not in rsID format - skipping SuSiE (PLINK requires rsIDs)")
        return(res_abf)
    }
    
    if (!exists("get_ld_matrix") || is.null(plink_bfile)) { 
        warning("LD tools missing. Skipping SuSiE."); 
        return(res_abf) 
    }
     
    ld_res <- get_ld_matrix(d1$snp, plink_bfile, plink_bin, keep_file = keep_file)
    
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

    # ==========================================================================
    # SuSiE-RSS Fine-Mapping: Set seeds for reproducibility
    # ==========================================================================
    # susie_rss() uses iterative variational inference which involves
    # stochastic initialization. We use offset seeds from the global seed
    # to ensure: (1) reproducibility across runs, (2) different seeds for
    # GWAS vs QTL to avoid artificial correlation.
    # Offsets: +1 for GWAS, +2 for QTL
    if (exists("global_seed", mode = "numeric")) {
        set.seed(global_seed + 1)
        message("[SuSiE] Seed set for GWAS fine-mapping (global_seed + 1)")
    }
    susie_1 <- try(susie_rss(bhat = d1_susie$beta, shat = sqrt(d1_susie$varbeta), R = R_susie, n = d1_susie$N[1], verbose = FALSE), silent = TRUE)

    if (exists("global_seed", mode = "numeric")) {
        set.seed(global_seed + 2)
        message("[SuSiE] Seed set for QTL fine-mapping (global_seed + 2)")
    }
    susie_2 <- try(susie_rss(bhat = d2_susie$beta, shat = sqrt(d2_susie$varbeta), R = R_susie, n = d2_susie$N[1], verbose = FALSE), silent = TRUE)
    
    if (inherits(susie_1, "try-error") || inherits(susie_2, "try-error")) { 
        warning("SuSiE failed to converge."); 
        return(res_abf) 
    }
    
    res_susie <- coloc.susie(susie_1, susie_2)
    message("  SuSiE Coloc complete.")
    
    res_abf$susie_result <- res_susie
    return(res_abf)
}
