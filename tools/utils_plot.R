# ------------------------------------------------------------------------------
# src/utils_plot.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(GenomicRanges)
  library(dplyr)
  library(data.table)
  library(vroom)
  library(glue)
  library(forcats)
  library(purrr)
  library(ggrepel)
  library(rtracklayer)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(tidyr)
})

# ==============================================================================
# Helper: Force Atomic (The "Nuclear Option" for Data Cleaning)
# ==============================================================================
resolve_col <- function(df, preferred, alternatives = NULL) {
  if (preferred %in% names(df)) return(preferred)
  if (!is.null(alternatives)) {
    for (alt in alternatives) {
      if (alt %in% names(df)) return(alt)
    }
  }
  available <- paste(head(names(df), 10), collapse=", ")
  stop(glue("Column '{preferred}' not found. Checked alternatives: {paste(alternatives, collapse=', ')}. Available cols: {available}..."))
}

# FORCE ATOMIC: Aggressively converts columns to basic types or drops them
force_atomic <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    
    # 1. Convert to base data.frame to strip tibble/grouped_df classes first
    df <- as.data.frame(df)
    
    # 2. Define safe types
    cols_to_keep <- c()
    
    for (col in names(df)) {
        val <- df[[col]]
        
        # Scenario A: It's already atomic (numeric, char, logical, integer) -> Keep
        if (is.vector(val) && is.atomic(val)) {
            cols_to_keep <- c(cols_to_keep, col)
            next
        }
        
        # Scenario B: It's an Rle (Run Length Encoding) -> Convert to vector
        if (inherits(val, "Rle")) {
            df[[col]] <- as.vector(val)
            # Re-check if atomic after conversion
            if (is.atomic(df[[col]])) {
                cols_to_keep <- c(cols_to_keep, col)
            }
            next
        }
        
        # Scenario C: List / SimpleList / AsIs
        if (is.list(val) || inherits(val, "List") || inherits(val, "SimpleList") || inherits(val, "AsIs")) {
            # Try to unlist if length matches
            try_flat <- tryCatch(unlist(val), error=function(e) NULL)
            
            if (!is.null(try_flat) && length(try_flat) == nrow(df)) {
                df[[col]] <- try_flat
                cols_to_keep <- c(cols_to_keep, col)
            } else {
                # Force to character string
                try_string <- tryCatch(vapply(val, function(x) paste(as.character(x), collapse=","), character(1)), error=function(e) NULL)
                if (!is.null(try_string) && length(try_string) == nrow(df)) {
                    df[[col]] <- try_string
                    cols_to_keep <- c(cols_to_keep, col)
                }
            }
        }
    }
    
    # 3. Return only safe columns
    return(df[, cols_to_keep, drop = FALSE])
}

# ==============================================================================
# 1. LD Calculation
# ==============================================================================

ld_extract <- function(variants, bfile, plink_bin) {
  shell <- ifelse(Sys.info()["sysname"] == "Mac", "cmd", "sh")
  fn <- tempfile()
  variants <- unique(as.character(variants))
  
  if (length(variants) == 0) return(NULL)

  data.frame(variants) %>%
    vroom::vroom_write(fn, col_names = FALSE)
  
  fun2 <- paste0(
    shQuote(plink_bin, type = shell),
    " --bfile ", shQuote(bfile, type = shell),
    " --extract ", shQuote(fn, type = shell),
    " --r inter-chr ",
    " --out ", shQuote(fn, type = shell)
  )
  
  system(fun2, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ld_file <- paste0(fn, ".ld")
  if (!file.exists(ld_file) || file.size(ld_file) == 0) {
      unlink(paste0(fn, "*"))
      return(NULL)
  }

  res <- tryCatch({
      data.table::fread(ld_file, header = TRUE)
  }, error = function(e) NULL)

  unlink(paste0(fn, "*"))

  if (is.null(res) || nrow(res) == 0) return(NULL)
  if (!all(c("SNP_A", "SNP_B", "R") %in% colnames(res))) return(NULL)

  res <- res %>% dplyr::select(SNP_A, SNP_B, R)
  res_inv <- res %>% 
      mutate(tmp_A = SNP_A, tmp_B = SNP_B) %>%
      mutate(SNP_A = tmp_B, SNP_B = tmp_A) %>%
      dplyr::select(SNP_A, SNP_B, R)

  return(bind_rows(res, res_inv) %>% distinct())
}

# ==============================================================================
# 2. LD Plot
# ==============================================================================

LD_plot <- function(df, ld_df = NULL, lead_snps = NULL, rsid = rsid, chromosome = chromosome, position = position, p_value = p_value,
                    p_value_threshold = 0.0000001, clump_kb = 1000, clump_r2 = 0.2, 
                    xlim = NULL,
                    bfile = NULL, plink_bin = NULL, plot_title = NULL, plot_subtitle = NULL, n_row = 2, region_recomb = region_recomb) {

  # CLEAN INPUT IMMEDIATELY using base R to avoid dplyr errors on dirty data
  df <- force_atomic(df)

  if (any(duplicated(colnames(df)))) {
      df <- df[, !duplicated(colnames(df))]
  }

  # Use base R subsetting to be safe before dplyr
  # Mapping: rsid, chromosome, position, p_value
  # We construct a clean dataframe
  clean_df <- data.frame(
      rsid = as.character(df[[resolve_col(df, "rsid", c("snp", "SNP", "marker", "id"))]]),
      chromosome = as.character(df[[resolve_col(df, "chromosome", c("chr", "CHR", "chrom"))]]),
      position = as.numeric(df[[resolve_col(df, "position", c("pos", "POS", "bp", "BP"))]]),
      p_value = as.numeric(df[[resolve_col(df, "p_value", c("pval", "P", "p", "P.gwas", "P.qtl"))]]),
      stringsAsFactors = FALSE
  )
  
  df <- clean_df

  if(nrow(df) == 0) return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No Data")))

  if (!is.null(xlim)) {
      df_view <- df %>% filter(position >= xlim[1] & position <= xlim[2])
  } else {
      df_view <- df
  }
  
  if (nrow(df_view) == 0) return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No SNPs in Window")))

  # Identify Lead SNPs
  if (!is.null(lead_snps)) {
    indep_snps <- df %>% dplyr::select(rsid, pval = p_value) %>% filter(rsid %in% lead_snps)
  } else {
    if (is.null(bfile) || is.null(plink_bin)) {
        indep_snps <- df %>% dplyr::select(rsid, pval = p_value) %>% arrange(pval) %>% slice(1)
    } else {
        indep_snps <- df %>% dplyr::select(rsid, pval = p_value) %>% filter(pval < p_value_threshold)
        if (nrow(indep_snps) == 0) {
            indep_snps <- df %>% arrange(pval) %>% slice(1) %>% dplyr::select(rsid, pval = p_value)
        }
    }
  }

  if(nrow(indep_snps) == 0) return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No Significant SNPs")))

  locus_snps <- df_view %>%
    dplyr::select(chromosome, position, rsid) %>%
    mutate(lead_rsid = indep_snps$rsid[1]) 

  # LD Calculation
  locus_snps_ld <- tryCatch({
      if (!is.null(ld_df)) {
          if (!all(c("SNP_A", "SNP_B") %in% colnames(ld_df))) stop("Invalid ld_df")
          if (!"r2" %in% colnames(ld_df) && "R" %in% colnames(ld_df)) ld_df <- ld_df %>% mutate(r2 = R^2)
          locus_snps %>% left_join(ld_df, by = c("lead_rsid" = "SNP_A", "rsid" = "SNP_B")) %>% dplyr::select(lead_rsid, rsid, r2)
      } else if (is.null(bfile) || is.null(plink_bin)) {
          locus_snps %>% mutate(r2 = NA) %>% dplyr::select(lead_rsid, rsid, r2)
      } else {
          snps_in_view <- unique(locus_snps$rsid)
          lead_snp <- indep_snps$rsid[1]
          if(is.na(lead_snp)) stop("No lead SNP")

          query_snps <- unique(c(lead_snp, snps_in_view))
          ld_res <- ld_extract(query_snps, bfile = bfile, plink_bin = plink_bin)
          
          if (is.null(ld_res) || nrow(ld_res) == 0) {
               locus_snps %>% mutate(r2 = NA) %>% dplyr::select(lead_rsid, rsid, r2)
          } else {
               ld_res %>% 
                 filter(SNP_A == lead_snp) %>%
                 mutate(r2 = abs(R)^2) %>%
                 dplyr::select(lead_rsid = SNP_A, rsid = SNP_B, r2)
          }
      }
  }, error = function(e) {
      locus_snps %>% mutate(r2 = NA) %>% dplyr::select(lead_rsid, rsid, r2)
  })

  RA_plot_data <- locus_snps %>% 
      left_join(locus_snps_ld, by=c("lead_rsid", "rsid")) %>%
      left_join(df_view %>% dplyr::select(rsid, p_value), by="rsid") %>%
      mutate(r2 = coalesce(as.numeric(r2), 0.1))

  RA_plot_data <- RA_plot_data %>% mutate(lead = (rsid == indep_snps$rsid[1])) %>%
    mutate(color_code = as.character(cut(r2, breaks = c(-0.1, 0.2, 0.4, 0.6, 0.8, 1.0), labels = c("blue4", "blue", "darkgreen", "orange", "red"), include.lowest = TRUE))) %>%
    mutate(color_code = ifelse(lead == TRUE, "purple", color_code))

  maxlogP <- max(-log10(RA_plot_data$p_value), na.rm=TRUE); if(is.infinite(maxlogP) || maxlogP==0) maxlogP <- 10

  if (!is.null(region_recomb) && nrow(region_recomb) > 0) {
      if(!is.null(xlim)) region_recomb <- region_recomb %>% filter(`Position(bp)` >= xlim[1] & `Position(bp)` <= xlim[2])
      recomb_max <- max(region_recomb$`Rate(cM/Mb)`, na.rm=TRUE)
      if(is.infinite(recomb_max) || recomb_max==0) recomb_max <- 100
      scale_factor <- maxlogP / recomb_max
      region_recomb_scaled <- region_recomb %>% mutate(scaled_rate = `Rate(cM/Mb)` * scale_factor)
  } else {
      scale_factor <- 1
      region_recomb_scaled <- data.frame(`Position(bp)` = numeric(0), scaled_rate = numeric(0))
  }

  p <- ggplot() +
    geom_point(data = RA_plot_data, mapping = aes(x = position/1e6, y = -log10(p_value), fill = color_code, size = lead, shape = lead), alpha = 0.8) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey") +
    scale_fill_identity(name = expression(LD ~ (r^2)),
                        labels = c("Lead SNP", "0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2"),
                        breaks = c("purple", "red", "orange", "darkgreen", "blue", "blue4"),
                        guide = guide_legend(override.aes = list(shape = 21, size = 5, color = "black", alpha = 1))) +
    scale_size_manual(values = c(3, 6), guide = "none") + scale_shape_manual(values = c(21, 23), guide = "none") +
    labs(title = plot_title, subtitle = plot_subtitle, x = "Position (Mb)") +
    theme_light(base_size = 14) + theme(legend.position = "right", plot.margin = margin(b=0, unit="cm"))

  if(!is.null(xlim)) {
      p <- p + scale_x_continuous(limits = xlim/1e6, n.breaks = 6, expand = c(0, 0))
  } else {
      p <- p + scale_x_continuous(n.breaks = 6)
  }

  if (nrow(region_recomb_scaled) > 0) {
     p <- p +
       geom_line(data = region_recomb_scaled, mapping = aes(x = `Position(bp)`/1e6, y = scaled_rate, linetype = "Recomb Rate"), color = "blue", linewidth = 0.8, alpha = 0.6) +
       scale_y_continuous(name = expression(-log[10]("p-value")), sec.axis = sec_axis(~ . / scale_factor, name = "Recomb Rate (cM/Mb)")) +
       scale_linetype_manual(name = "", values = c("Recomb Rate" = 1), guide = "legend")
  } else {
     p <- p + scale_y_continuous(name = expression(-log[10]("p-value")))
  }

  return(p)
}

# ==============================================================================
# 2. Gene Track (Synced, Filtered & Safe - Base R Processing)
# ==============================================================================
genetrack <- function(chrom_str, xlim_bp, gtf_path = NULL) {
    chrom_clean <- gsub("^chr", "", as.character(chrom_str))
    chrom_str <- paste0("chr", chrom_clean)
    
    query_start <- xlim_bp[1]
    query_end <- xlim_bp[2]
    gr <- GRanges(seqnames = chrom_str, ranges = IRanges(start = query_start, end = query_end))
    
    tx_df <- data.frame(); exon_plot_data <- data.frame()
    
    if (!is.null(gtf_path) && file.exists(gtf_path)) {
        if (!requireNamespace("rtracklayer", quietly = TRUE)) return(ggplot() + theme_void() + geom_text(aes(x=0, y=0.5, label="rtracklayer missing")))
        tryCatch({
            subset_gtf <- rtracklayer::import(gtf_path, which = gr, feature.type = "exon")
            if (length(subset_gtf) > 0) {
                # Convert & Clean IMMEDIATELY using force_atomic
                df_gtf <- force_atomic(as.data.frame(subset_gtf))
                
                # Check critical columns
                if (!"transcript_id" %in% colnames(df_gtf)) df_gtf$transcript_id <- df_gtf$gene_id
                if ("gene_name" %in% colnames(df_gtf)) df_gtf$display_name <- df_gtf$gene_name else df_gtf$display_name <- df_gtf$gene_id
                
                # Filter Coding using base R subsetting (Safe)
                if ("gene_type" %in% colnames(df_gtf)) {
                     df_gtf <- df_gtf[df_gtf$gene_type == "protein_coding", , drop=FALSE]
                } else if ("gene_biotype" %in% colnames(df_gtf)) {
                     df_gtf <- df_gtf[df_gtf$gene_biotype == "protein_coding", , drop=FALSE]
                }
                
                if (nrow(df_gtf) == 0) return(ggplot() + theme_void() + geom_text(aes(x=mean(xlim_bp), y=0.5, label="No coding genes")))

                # Filter Longest Tx using Base R to avoid slice error
                # 1. Calc lengths
                df_gtf$len <- df_gtf$end - df_gtf$start
                # 2. Aggregate
                agg_res <- aggregate(len ~ transcript_id + display_name, data = df_gtf, FUN = sum)
                # 3. Find max per gene
                best_txs <- do.call(rbind, lapply(split(agg_res, agg_res$display_name), function(x) {
                    x[which.max(x$len), ]
                }))
                
                df_gtf_filtered <- df_gtf[df_gtf$transcript_id %in% best_txs$transcript_id, ]
                
                if(nrow(df_gtf_filtered) > 0) {
                    # Manual aggregation for tx_df
                    tx_df <- aggregate(cbind(start, end) ~ transcript_id + display_name + strand, data = df_gtf_filtered, 
                                       FUN = function(x) c(min(x), max(x)))
                    # Fix structure after aggregate with matrix output
                    tx_df <- data.frame(
                        tx_id = tx_df$transcript_id,
                        symbol = tx_df$display_name,
                        strand = tx_df$strand,
                        start = tx_df$start[,1], # min
                        end = tx_df$end[,2] # max
                    )
                    
                    exon_plot_data <- df_gtf_filtered[, c("start", "end", "transcript_id", "strand")]
                    colnames(exon_plot_data)[3] <- "tx_id"
                }
            }
        }, error = function(e) warning(paste("GTF error:", e$message)))
    } else {
        # TxDb Logic
        if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="TxDb missing")))
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
        subset_tx <- suppressWarnings(GenomicFeatures::transcriptsByOverlaps(txdb, gr))
        
        if (length(subset_tx) > 0) {
            symbols_vec <- tryCatch({ suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys=subset_tx$gene_id, column="SYMBOL", keytype="ENTREZID", multiVals="first")) }, error=function(e) NULL)
            
            subset_tx$symbol <- if(!is.null(symbols_vec)) symbols_vec else subset_tx$tx_name
            subset_tx$symbol[is.na(subset_tx$symbol)] <- subset_tx$tx_name[is.na(subset_tx$symbol)]
            
            # Convert & Clean
            subset_tx_df <- force_atomic(as.data.frame(subset_tx))
            subset_tx_df$width <- subset_tx_df$end - subset_tx_df$start
            
            # Heuristic Coding Filter (Base R)
            if ("symbol" %in% colnames(subset_tx_df)) {
                keep_idx <- !is.na(subset_tx_df$symbol) & 
                            !grepl("^LOC|^LINC|^MIR|^SNOR|^ENSG|^RN7|^RNU", subset_tx_df$symbol) &
                            subset_tx_df$symbol != ""
                subset_tx_df <- subset_tx_df[keep_idx, , drop=FALSE]
            }
            
            if (nrow(subset_tx_df) == 0) return(ggplot() + theme_void() + geom_text(aes(x=mean(xlim_bp), y=0.5, label="No coding genes")))

            # Find longest isoform (Base R)
            best_txs <- do.call(rbind, lapply(split(subset_tx_df, subset_tx_df$symbol), function(x) {
                x[which.max(x$width), ]
            }))
            best_tx_ids <- best_txs$tx_id
            
            if (length(best_tx_ids) > 0) {
                exons_all <- exons(txdb, filter=list(tx_id = best_tx_ids), columns=c("tx_id"))
                exon_plot_data <- force_atomic(as.data.frame(exons_all))
                
                # Manual unnest
                if (is.list(exon_plot_data$tx_id)) {
                    # Very primitive unnest to avoid dplyr dependency here
                    lens <- vapply(exon_plot_data$tx_id, length, integer(1))
                    ids_flat <- unlist(exon_plot_data$tx_id)
                    exon_plot_data <- exon_plot_data[rep(1:nrow(exon_plot_data), lens), ]
                    exon_plot_data$tx_id <- as.character(ids_flat)
                } else {
                    exon_plot_data$tx_id <- as.character(exon_plot_data$tx_id)
                }
                
                tx_df <- subset_tx_df[subset_tx_df$tx_id %in% best_tx_ids, ]
                tx_df <- tx_df[, c("tx_id", "start", "end", "symbol", "strand")]
                
                # Join strand to exons
                tx_lookup <- tx_df[, c("tx_id", "strand")]
                exon_plot_data <- merge(exon_plot_data, tx_lookup, by="tx_id")
            }
        }
    }
    
    if (nrow(tx_df) == 0) return(ggplot() + theme_void() + geom_text(aes(x=mean(xlim_bp), y=0.5, label="No coding genes")))
    
    # Clean up again just in case
    tx_df <- force_atomic(tx_df)
    exon_plot_data <- force_atomic(exon_plot_data)
    
    tx_df <- tx_df[order(tx_df$start), ]
    
    # Layout Logic
    levels <- rep(0, nrow(tx_df))
    if(nrow(tx_df) > 0) {
      ends <- c()
      for (i in 1:nrow(tx_df)) {
        found_level <- FALSE
        for (lvl in 1:length(ends)) {
          # Check for overlap + margin
          if (tx_df$start[i] > ends[lvl] + (xlim_bp[2]-xlim_bp[1])*0.05) { 
             levels[i] <- lvl
             ends[lvl] <- tx_df$end[i]
             found_level <- TRUE
             break
          }
        }
        if (!found_level) {
          levels[i] <- length(ends) + 1
          ends <- c(ends, tx_df$end[i])
        }
      }
    }
    tx_df$y <- levels
    
    exon_plot_data <- merge(exon_plot_data, tx_df[, c("tx_id", "y")], by="tx_id")
    strand_colors <- c("+" = "#D55E00", "-" = "#0072B2", "*" = "grey50")
    tx_df$label_text <- paste0(tx_df$symbol, ifelse(tx_df$strand == "+", " >", " <"))

    ggplot() + 
       geom_segment(data=tx_df, aes(x=start, xend=end, y=y, yend=y, color=strand), linewidth=0.5) +
       geom_rect(data=exon_plot_data, aes(xmin=start, xmax=end, ymin=y-0.25, ymax=y+0.25, fill=strand), color=NA) +
       geom_text(data=tx_df, aes(x=(start+end)/2, y=y, label=label_text), 
                       size=2.5, fontface="italic", vjust = -1.2) +
       
       scale_x_continuous(limits = xlim_bp, expand = c(0, 0)) + 
       scale_y_continuous(limits = c(0.5, max(tx_df$y) + 1), expand = c(0, 0)) + 
       scale_color_manual(values = strand_colors, guide="none") +
       scale_fill_manual(values = strand_colors, guide="none") +
       theme_void() + 
       theme(plot.margin = margin(t=0.1, b=0, unit="cm"))
}

# ==============================================================================
# 3. Composite Functions
# ==============================================================================

plot_qtl_association <- function(qtl_all_chrom, qtl_all_pvalue, leadSNP_DF, ld_df = NULL, gtf_path = NULL, region_recomb = NULL) {
    # CLEAN INPUT: Convert to clean DF first
    leadSNP_DF <- force_atomic(as.data.frame(leadSNP_DF))
    leadSNP_DF <- leadSNP_DF[, !duplicated(colnames(leadSNP_DF))]

    target_bp <- resolve_col(leadSNP_DF, "POS.qtl", c("POS", "pos", "BP", "bp", "POS.gwas"))
    target_chrom <- resolve_col(leadSNP_DF, qtl_all_chrom, c("CHR.qtl", "chr", "CHR", "chrom"))
    target_p <- resolve_col(leadSNP_DF, qtl_all_pvalue, c("P.qtl", "pval", "p", "P"))
    target_snp <- resolve_col(leadSNP_DF, "snp", c("SNP", "rsid", "RSID", "variant_id", "marker"))

    eQTL_leadSNP_DF <- leadSNP_DF %>%
        dplyr::select(rsid = all_of(target_snp),
                      chromosome = all_of(target_chrom),
                      position = all_of(target_bp),
                      p_value = all_of(target_p)) %>%
        force_atomic() 

    lead_snp_val <- if(exists("lead_SNP", envir=.GlobalEnv)) get("lead_SNP", envir=.GlobalEnv) else "LeadSNP"
    gene_sym     <- if(exists("geneSymbol", envir=.GlobalEnv)) get("geneSymbol", envir=.GlobalEnv) else "Gene"
    tissue_val   <- if(exists("type", envir=.GlobalEnv)) get("type", envir=.GlobalEnv) else "Type"
    plink_bfile_val <- if(exists("plink_bfile", envir=.GlobalEnv)) get("plink_bfile", envir=.GlobalEnv) else NULL

    # DEFINE VIEW WINDOW: +/- 200kb
    # Base R sort to avoid slice error
    eQTL_leadSNP_DF <- eQTL_leadSNP_DF[order(eQTL_leadSNP_DF$p_value), ]
    lead_snp_pos <- eQTL_leadSNP_DF$position[1]
    
    radius <- 200000 
    xlim_window <- c(lead_snp_pos - radius, lead_snp_pos + radius)

    p_assoc <- tryCatch({
        LD_plot(
            df = eQTL_leadSNP_DF, ld_df = ld_df, lead_snps = lead_snp_val, bfile = plink_bfile_val, plink_bin = "plink",
            xlim = xlim_window,
            plot_title = paste(lead_snp_val, gene_sym, "QTL"), plot_subtitle = "QTL Association (hg38)", region_recomb = region_recomb
        )
    }, error = function(e) {
        ggplot() + theme_void() + geom_text(aes(x=0, y=0, label=paste("Error:", e$message)))
    })

    p_track <- genetrack(eQTL_leadSNP_DF$chromosome[1], xlim_bp = xlim_window, gtf_path = gtf_path)
    ggarrange(p_assoc, p_track, ncol = 1, heights = c(3, 1), align = "v")
}

plot_gwas_association <- function(colocInputFile, qtl_all_chrom, leadSNP_DF, ld_df = NULL, gtf_path = NULL, region_recomb = NULL) {
    # CLEAN INPUT
    leadSNP_DF <- force_atomic(as.data.frame(leadSNP_DF))
    leadSNP_DF <- leadSNP_DF[, !duplicated(colnames(leadSNP_DF))]

    target_bp <- resolve_col(leadSNP_DF, "POS.qtl", c("POS.gwas", "POS", "pos", "BP"))
    target_chrom <- resolve_col(leadSNP_DF, qtl_all_chrom, c("CHR.gwas", "CHR", "chr"))
    target_p  <- resolve_col(leadSNP_DF, "P.gwas", c("P", "p", "pval"))
    target_snp <- resolve_col(leadSNP_DF, "snp", c("SNP", "rsid", "RSID", "variant_id", "marker"))

    trait_leadSNP_DF <- leadSNP_DF %>%
        dplyr::select(rsid = all_of(target_snp),
                      chromosome = all_of(target_chrom),
                      position = all_of(target_bp),
                      p_value = all_of(target_p)) %>%
        force_atomic()

    lead_snp_val <- if(exists("lead_SNP", envir=.GlobalEnv)) get("lead_SNP", envir=.GlobalEnv) else "LeadSNP"
    trait_name   <- if(exists("trait", envir=.GlobalEnv)) get("trait", envir=.GlobalEnv) else "Trait"
    plink_bfile_val <- if(exists("plink_bfile", envir=.GlobalEnv)) get("plink_bfile", envir=.GlobalEnv) else NULL
    
    # Base R sort
    trait_leadSNP_DF <- trait_leadSNP_DF[order(trait_leadSNP_DF$p_value), ]
    lead_snp_pos <- trait_leadSNP_DF$position[1]
    
    radius <- 200000 
    xlim_window <- c(lead_snp_pos - radius, lead_snp_pos + radius)

    p_assoc <- tryCatch({
        LD_plot(df = trait_leadSNP_DF, ld_df = ld_df, lead_snps = lead_snp_val, bfile = plink_bfile_val, plink_bin = "plink", 
                xlim = xlim_window,
                plot_title = paste(lead_snp_val, trait_name, "GWAS"), region_recomb = region_recomb)
    }, error = function(e) {
        ggplot() + theme_void() + geom_text(aes(x=0, y=0, label=paste("Error:", e$message)))
    })

    p_track <- genetrack(trait_leadSNP_DF$chromosome[1], xlim_bp = xlim_window, gtf_path = gtf_path)
    ggarrange(p_assoc, p_track, ncol = 1, heights = c(3, 1), align = "v")
}
