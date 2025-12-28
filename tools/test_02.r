# ------------------------------------------------------------------------------
# utils_plot.R - Final Stable Version
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
  library(grid)
})

# ==============================================================================
# Helper Functions
# ==============================================================================

resolve_col <- function(df, preferred, alternatives = NULL) {
  if (preferred %in% names(df)) return(preferred)
  if (!is.null(alternatives)) {
    for (alt in alternatives) {
      if (alt %in% names(df)) return(alt)
    }
  }
  # Case-insensitive fallback
  all_cols <- names(df)
  match_idx <- which(tolower(all_cols) == tolower(preferred))
  if (length(match_idx) > 0) return(all_cols[match_idx[1]])
  
  if (!is.null(alternatives)) {
    for (alt in alternatives) {
      match_idx <- which(tolower(all_cols) == tolower(alt))
      if (length(match_idx) > 0) return(all_cols[match_idx[1]])
    }
  }
  return(NULL)
}

safe_df_convert <- function(df) {
    if (is.null(df)) return(NULL)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (nrow(df) == 0) return(NULL)
    return(df)
}

safe_numeric <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0) return(default)
    suppressWarnings(as.numeric(as.character(x)))
}

safe_character <- function(x, default = NA_character_) {
    if (is.null(x) || length(x) == 0) return(default)
    as.character(x)
}

# ==============================================================================
# LD Calculation
# ==============================================================================

ld_extract <- function(variants, bfile, plink_bin) {
  variants <- unique(as.character(variants))
  if (length(variants) == 0) return(NULL)

  fn <- tempfile()
  on.exit(unlink(paste0(fn, "*")), add = TRUE)

  tryCatch({
    data.frame(variants, stringsAsFactors = FALSE) %>%
      vroom::vroom_write(fn, col_names = FALSE)

    shell <- ifelse(Sys.info()["sysname"] == "Darwin", "cmd", "sh")
    cmd <- paste0(
      shQuote(plink_bin, type = shell),
      " --bfile ", shQuote(bfile, type = shell),
      " --extract ", shQuote(fn, type = shell),
      " --r inter-chr --out ", shQuote(fn, type = shell)
    )

    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    ld_file <- paste0(fn, ".ld")
    if (!file.exists(ld_file) || file.size(ld_file) == 0) return(NULL)

    res <- data.table::fread(ld_file, header = TRUE)
    if (!all(c("SNP_A", "SNP_B", "R") %in% colnames(res))) return(NULL)

    res <- res %>% dplyr::select(SNP_A, SNP_B, R)
    res_inv <- data.frame(SNP_A = res$SNP_B, SNP_B = res$SNP_A, R = res$R)

    return(bind_rows(res, res_inv) %>% distinct())
  }, error = function(e) {
    message(glue("[LD] Extraction failed: {e$message}"))
    return(NULL)
  })
}

# ==============================================================================
# LD Plot (Fixed)
# ==============================================================================

LD_plot <- function(df, ld_df = NULL, lead_snps = NULL,
                    bfile = NULL, plink_bin = "plink",
                    plot_title = NULL, plot_subtitle = NULL,
                    region_recomb = NULL, xlim = NULL) {

  df <- safe_df_convert(df)
  if (is.null(df)) {
    return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No Input Data")))
  }

  col_snp <- resolve_col(df, "rsid", c("snp", "SNP", "SNPID", "marker"))
  col_chr <- resolve_col(df, "chromosome", c("CHR", "chr", "CHR.qtl", "CHR.gwas"))
  col_pos <- resolve_col(df, "position", c("POS", "pos", "POS.qtl", "POS.gwas", "BP"))
  col_p   <- resolve_col(df, "p_value", c("P", "pval", "P.qtl", "P.gwas"))

  if (is.null(col_snp) || is.null(col_chr) || is.null(col_pos) || is.null(col_p)) {
    return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="Missing Columns")))
  }

  plot_df <- data.frame(
    rsid = safe_character(df[[col_snp]]),
    chromosome = safe_character(df[[col_chr]]),
    position = safe_numeric(df[[col_pos]]),
    p_value = safe_numeric(df[[col_p]]),
    stringsAsFactors = FALSE
  )

  plot_df <- plot_df[
    !is.na(plot_df$rsid) & plot_df$rsid != "" &
    !is.na(plot_df$position) &
    !is.na(plot_df$p_value) & plot_df$p_value > 0,
  ]

  if (nrow(plot_df) == 0) {
    return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No Valid SNPs")))
  }

  # Apply xlim filter if provided
  if (!is.null(xlim) && length(xlim) == 2) {
    plot_df <- plot_df[plot_df$position >= xlim[1] & plot_df$position <= xlim[2], ]
    if (nrow(plot_df) == 0) {
      return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No SNPs in Window")))
    }
  }

  message(glue("[LD_plot] Plotting {nrow(plot_df)} SNPs"))

  # Determine lead SNP
  if (!is.null(lead_snps) && length(lead_snps) > 0) {
    lead_snp <- lead_snps[1]
    if (!lead_snp %in% plot_df$rsid) {
      lead_snp <- plot_df$rsid[which.min(plot_df$p_value)]
    }
  } else {
    lead_snp <- plot_df$rsid[which.min(plot_df$p_value)]
  }

  message(glue("[LD_plot] Lead SNP: {lead_snp}"))

  # Initialize r2 column
  plot_df$r2 <- NA_real_

  # Calculate LD if bfile exists
  if (!is.null(bfile) && file.exists(paste0(bfile, ".bed"))) {
    ld_result <- ld_extract(unique(c(lead_snp, plot_df$rsid)), bfile, plink_bin)

    if (!is.null(ld_result)) {
      ld_to_lead <- ld_result[ld_result$SNP_A == lead_snp, ]

      if (nrow(ld_to_lead) > 0) {
        ld_to_lead$r2_calc <- abs(ld_to_lead$R)^2

        # Merge using specific temp column to avoid conflicts
        plot_df <- merge(
          plot_df,
          ld_to_lead[, c("SNP_B", "r2_calc")],
          by.x = "rsid",
          by.y = "SNP_B",
          all.x = TRUE
        )

        # Assign calculated LD, keeping NA for unmatched
        plot_df$r2 <- ifelse(!is.na(plot_df$r2_calc), plot_df$r2_calc, plot_df$r2)
        plot_df$r2_calc <- NULL
      }
    }
  }

  # Lead SNP always r2=1, background=0.1
  plot_df$r2[plot_df$rsid == lead_snp] <- 1.0
  plot_df$r2[is.na(plot_df$r2)] <- 0.1

  high_ld_count <- sum(plot_df$r2 > 0.2 & plot_df$rsid != lead_snp)
  message(glue("[LD_plot] Found {high_ld_count} proxy SNPs with r2 > 0.2"))

  plot_df$is_lead <- plot_df$rsid == lead_snp

  plot_df$color_code <- cut(
    plot_df$r2,
    breaks = c(-0.1, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c("navy", "skyblue", "springgreen3", "orange", "red"),
    include.lowest = TRUE
  )
  plot_df$color_code <- as.character(plot_df$color_code)
  plot_df$color_code[plot_df$is_lead] <- "purple"

  maxlogP <- max(-log10(plot_df$p_value), na.rm = TRUE)
  if (is.infinite(maxlogP) || maxlogP < 1) maxlogP <- 10

  p <- ggplot(plot_df, aes(x = position/1e6, y = -log10(p_value))) +
    geom_point(aes(fill = color_code, size = is_lead, shape = is_lead), 
               alpha = 0.85, stroke = 0.3, color = "black") +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey50") +
    scale_fill_identity(
      name = expression(LD ~ (r^2)),
      labels = c("Lead SNP", "0.8-1.0", "0.6-0.8", "0.4-0.6", "0.2-0.4", "0-0.2"),
      breaks = c("purple", "red", "orange", "springgreen3", "skyblue", "navy"),
      guide = guide_legend(override.aes = list(shape = 21, size = 4, alpha = 1))
    ) +
    scale_size_manual(values = c(2.5, 5.5), guide = "none") +
    scale_shape_manual(values = c(21, 23), guide = "none") +
    labs(title = plot_title, subtitle = plot_subtitle, 
         x = NULL, y = expression(-log[10]("p-value"))) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(b = 0, unit = "pt"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  if (!is.null(xlim) && length(xlim) == 2) {
    p <- p + scale_x_continuous(limits = xlim/1e6, expand = c(0, 0))
  }

  # Add recombination rate
  if (!is.null(region_recomb) && is.data.frame(region_recomb) && nrow(region_recomb) > 0) {
    r_pos_col <- resolve_col(region_recomb, "Position(bp)", c("position", "pos", "bp"))
    r_rate_col <- resolve_col(region_recomb, "Rate(cM/Mb)", c("rate", "recomb", "cM/Mb"))

    if (!is.null(r_pos_col) && !is.null(r_rate_col)) {
      region_recomb$pos_mapped <- safe_numeric(region_recomb[[r_pos_col]])
      region_recomb$rate_mapped <- safe_numeric(region_recomb[[r_rate_col]])
      
      region_recomb <- region_recomb[!is.na(region_recomb$pos_mapped) & !is.na(region_recomb$rate_mapped), ]
      
      if (!is.null(xlim) && length(xlim) == 2) {
        region_recomb <- region_recomb[
          region_recomb$pos_mapped >= xlim[1] & 
          region_recomb$pos_mapped <= xlim[2], 
        ]
      }

      if (nrow(region_recomb) > 0) {
        recomb_max <- max(region_recomb$rate_mapped, na.rm = TRUE)
        if (!is.infinite(recomb_max) && recomb_max > 0) {
          scale_factor <- (maxlogP * 0.9) / recomb_max
          region_recomb$scaled_rate <- region_recomb$rate_mapped * scale_factor

          p <- p +
            geom_area(data = region_recomb,
                      aes(x = pos_mapped/1e6, y = scaled_rate),
                      fill = "#E6F5FF", alpha = 0.6, inherit.aes = FALSE) +
            geom_line(data = region_recomb,
                      aes(x = pos_mapped/1e6, y = scaled_rate),
                      color = "#4682B4", linewidth = 0.7, inherit.aes = FALSE) +
            scale_y_continuous(
              name = expression(-log[10]("p-value")),
              sec.axis = sec_axis(~ . / scale_factor, name = "Recomb Rate (cM/Mb)")
            )
        }
      }
    }
  }

  return(p)
}

# ==============================================================================
# Gene Track (Simplified & Stable)
# ==============================================================================

genetrack <- function(chrom_str, xlim_bp, gtf_path = NULL) {
    chrom_clean <- gsub("^chr", "", as.character(chrom_str))
    chrom_str <- paste0("chr", chrom_clean)
    
    if (is.null(xlim_bp) || length(xlim_bp) != 2) {
      xlim_bp <- c(0, 1e6)
    }

    BPStart <- xlim_bp[1]
    BPStop <- xlim_bp[2]

    gr <- GRanges(seqnames = chrom_str, ranges = IRanges(start = BPStart, end = BPStop))
    tx_df <- NULL

    # Try GTF first
    if (!is.null(gtf_path) && file.exists(gtf_path)) {
        tryCatch({
            subset_gtf <- rtracklayer::import(gtf_path, which = gr, feature.type = "exon")
            if (length(subset_gtf) > 0) {
                df_gtf <- as.data.frame(subset_gtf)

                if ("gene_name" %in% names(df_gtf)) {
                    df_gtf$symbol <- df_gtf$gene_name
                } else if ("gene_id" %in% names(df_gtf)) {
                    df_gtf$symbol <- df_gtf$gene_id
                } else {
                    df_gtf$symbol <- "Gene"
                }

                if ("gene_type" %in% names(df_gtf)) {
                    df_gtf <- df_gtf[df_gtf$gene_type == "protein_coding", ]
                }

                if (nrow(df_gtf) > 0) {
                    tx_summary <- df_gtf %>%
                        group_by(symbol) %>%
                        summarize(
                            start = min(start),
                            end = max(end),
                            .groups = 'drop'
                        )
                    tx_df <- as.data.frame(tx_summary)
                }
            }
        }, error = function(e) {
            message(glue("[Gene Track] GTF error: {e$message}"))
        })
    }

    # Fallback to TxDb
    if (is.null(tx_df) || nrow(tx_df) == 0) {
        tryCatch({
            txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
            subset_tx <- GenomicFeatures::transcriptsByOverlaps(txdb, gr)

            if (length(subset_tx) > 0) {
                symbols <- suppressMessages(
                    AnnotationDbi::mapIds(
                        org.Hs.eg.db::org.Hs.eg.db,
                        keys = subset_tx$gene_id,
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first"
                    )
                )

                tx_df <- data.frame(
                    symbol = ifelse(is.na(symbols), subset_tx$tx_name, symbols),
                    start = start(subset_tx),
                    end = end(subset_tx),
                    stringsAsFactors = FALSE
                )

                tx_df <- tx_df[!grepl("^(LOC|LINC|MIR|SNOR|RN7|RNU)", tx_df$symbol), ]

                if (nrow(tx_df) > 0) {
                    tx_df <- tx_df %>%
                        group_by(symbol) %>%
                        summarize(start = min(start), end = max(end), .groups = 'drop') %>%
                        as.data.frame()
                }
            }
        }, error = function(e) {
            message(glue("[Gene Track] TxDb error: {e$message}"))
        })
    }

    if (is.null(tx_df) || nrow(tx_df) == 0) {
        return(
          ggplot() + theme_void() +
          geom_text(aes(x = mean(xlim_bp)/1e6, y = 0.5, label = "No genes"), 
                    size = 4, color = "grey50") +
          scale_x_continuous(limits = xlim_bp/1e6, expand = c(0, 0)) +
          theme(
            plot.margin = margin(t = 5, b = 10, unit = "pt"),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.ticks.x = element_line(color = "black", linewidth = 0.5),
            axis.line.x = element_line(color = "black", linewidth = 0.5)
          ) +
          labs(x = paste0("Position on ", chrom_str, " (Mb)"))
        )
    }

    message(glue("[Gene Track] Drawing {nrow(tx_df)} genes"))

    tx_df <- tx_df[order(tx_df$start), ]

    # Simple staggered layout
    tx_df$y <- 1
    if (nrow(tx_df) > 1) {
        for (i in 2:nrow(tx_df)) {
            overlap_margin <- (BPStop - BPStart) * 0.02
            if (tx_df$start[i] < (tx_df$end[i-1] + overlap_margin)) {
                tx_df$y[i] <- ifelse(tx_df$y[i-1] == 1, 2, 1)
            }
        }
    }

    ggplot(tx_df) +
        geom_segment(aes(x = start/1e6, xend = end/1e6, y = y, yend = y),
                     color = "steelblue", linewidth = 2) +
        geom_text(aes(x = (start + end)/2e6, y = y, label = symbol),
                  size = 3.2, fontface = "bold.italic", vjust = -0.8) +
        scale_x_continuous(limits = xlim_bp/1e6, expand = c(0, 0)) +
        scale_y_continuous(limits = c(0.5, 2.5), expand = c(0, 0)) +
        theme_void() +
        theme(
          plot.margin = margin(t = 5, b = 10, unit = "pt"),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.ticks.x = element_line(color = "black", linewidth = 0.5),
          axis.line.x = element_line(color = "black", linewidth = 0.5)
        ) +
        labs(x = paste0("Position on ", chrom_str, " (Mb)"))
}

# ==============================================================================
# Composite Functions
# ==============================================================================

plot_qtl_association <- function(qtl_all_chrom, qtl_all_pvalue, leadSNP_DF,
                                 ld_df = NULL, gtf_path = NULL, region_recomb = NULL) {

    message("[plot_qtl_association] Starting...")

    leadSNP_DF <- safe_df_convert(leadSNP_DF)
    if (is.null(leadSNP_DF)) return(NULL)

    lead_snp_val <- if(exists("lead_SNP", envir=.GlobalEnv)) get("lead_SNP", envir=.GlobalEnv) else NULL
    gene_sym <- if(exists("geneSymbol", envir=.GlobalEnv)) get("geneSymbol", envir=.GlobalEnv) else "Gene"
    plink_bfile_val <- if(exists("plink_bfile", envir=.GlobalEnv)) get("plink_bfile", envir=.GlobalEnv) else NULL

    pos_col <- resolve_col(leadSNP_DF, "POS.qtl", c("POS", "pos", "BP", "POS.gwas", "position"))
    chr_col <- resolve_col(leadSNP_DF, qtl_all_chrom, c("CHR", "CHR.qtl", "chr", "chromosome"))

    if (is.null(pos_col) || is.null(chr_col)) {
      message("[ERROR] Cannot resolve position/chr columns")
      return(NULL)
    }

    positions <- safe_numeric(leadSNP_DF[[pos_col]])
    positions <- positions[!is.na(positions) & positions > 0]

    if (!is.null(lead_snp_val)) {
        snp_col <- resolve_col(leadSNP_DF, "rsid", c("snp", "SNP", "marker"))
        if (!is.null(snp_col)) {
          lead_row <- leadSNP_DF[leadSNP_DF[[snp_col]] == lead_snp_val, ]
          if (nrow(lead_row) > 0) {
              center_pos <- safe_numeric(lead_row[[pos_col]][1])
              min_pos <- center_pos - 200000
              max_pos <- center_pos + 200000
          } else {
              min_pos <- min(positions)
              max_pos <- max(positions)
          }
        } else {
          min_pos <- min(positions)
          max_pos <- max(positions)
        }
    } else {
        min_pos <- min(positions)
        max_pos <- max(positions)
    }

    chrom_num <- safe_character(leadSNP_DF[[chr_col]])[1]

    p_assoc <- tryCatch({
        LD_plot(
            df = leadSNP_DF,
            ld_df = ld_df,
            lead_snps = lead_snp_val,
            bfile = plink_bfile_val,
            plink_bin = "plink",
            plot_title = paste(lead_snp_val, gene_sym, "QTL"),
            plot_subtitle = "QTL Association (hg38)",
            region_recomb = region_recomb,
            xlim = c(min_pos, max_pos)
        )
    }, error = function(e) {
        message(glue("[ERROR] LD plot failed: {e$message}"))
        ggplot() + theme_void() + geom_text(aes(x=0, y=0, label=paste("Error:", e$message)))
    })

    p_track <- tryCatch({
        genetrack(chrom_num, c(min_pos, max_pos), gtf_path = gtf_path)
    }, error = function(e) {
        message(glue("[ERROR] Gene track failed: {e$message}"))
        ggplot() + theme_void()
    })

    ggarrange(p_assoc, p_track, ncol = 1, heights = c(2, 1), align = "v")
}

plot_gwas_association <- function(colocInputFile, qtl_all_chrom, leadSNP_DF,
                                  ld_df = NULL, gtf_path = NULL, region_recomb = NULL) {

    message("[plot_gwas_association] Starting...")

    leadSNP_DF <- safe_df_convert(leadSNP_DF)
    if (is.null(leadSNP_DF)) return(NULL)

    lead_snp_val <- if(exists("lead_SNP", envir=.GlobalEnv)) get("lead_SNP", envir=.GlobalEnv) else NULL
    trait_name <- if(exists("trait", envir=.GlobalEnv)) get("trait", envir=.GlobalEnv) else "GWAS"
    plink_bfile_val <- if(exists("plink_bfile", envir=.GlobalEnv)) get("plink_bfile", envir=.GlobalEnv) else NULL

    pos_col <- resolve_col(leadSNP_DF, "POS.gwas", c("POS", "pos", "BP", "POS.qtl", "position"))
    chr_col <- resolve_col(leadSNP_DF, "CHR.gwas", c("CHR", "chr", "CHR.qtl", "chromosome"))

    if (is.null(pos_col) || is.null(chr_col)) {
      message("[ERROR] Cannot resolve position/chr columns")
      return(NULL)
    }

    positions <- safe_numeric(leadSNP_DF[[pos_col]])
    positions <- positions[!is.na(positions) & positions > 0]

    if (!is.null(lead_snp_val)) {
        snp_col <- resolve_col(leadSNP_DF, "rsid", c("snp", "SNP", "marker"))
        if (!is.null(snp_col)) {
          lead_row <- leadSNP_DF[leadSNP_DF[[snp_col]] == lead_snp_val, ]
          if (nrow(lead_row) > 0) {
              center_pos <- safe_numeric(lead_row[[pos_col]][1])
              min_pos <- center_pos - 200000
              max_pos <- center_pos + 200000
          } else {
              min_pos <- min(positions)
              max_pos <- max(positions)
          }
        } else {
          min_pos <- min(positions)
          max_pos <- max(positions)
        }
    } else {
        min_pos <- min(positions)
        max_pos <- max(positions)
    }

    chrom_num <- safe_character(leadSNP_DF[[chr_col]])[1]

    p_assoc <- tryCatch({
        LD_plot(
            df = leadSNP_DF,
            ld_df = ld_df,
            lead_snps = lead_snp_val,
            bfile = plink_bfile_val,
            plink_bin = "plink",
            plot_title = paste(lead_snp_val, trait_name, "GWAS"),
            plot_subtitle = "GWAS Association (hg38)",
            region_recomb = region_recomb,
            xlim = c(min_pos, max_pos)
        )
    }, error = function(e) {
        message(glue("[ERROR] LD plot failed: {e$message}"))
        ggplot() + theme_void() + geom_text(aes(x=0, y=0, label=paste("Error:", e$message)))
    })

    p_track <- tryCatch({
        genetrack(chrom_num, c(min_pos, max_pos), gtf_path = gtf_path)
    }, error = function(e) {
        message(glue("[ERROR] Gene track failed: {e$message}"))
        ggplot() + theme_void()
    })

    ggarrange(p_assoc, p_track, ncol = 1, heights = c(2, 1), align = "v")
}
