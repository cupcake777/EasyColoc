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
  library(grid)
})

resolve_col <- function(df, preferred, alternatives = NULL) {
  if (preferred %in% names(df)) return(preferred)
  if (!is.null(alternatives)) {
    for (alt in alternatives) {
      if (alt %in% names(df)) return(alt)
    }
  }
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

ld_extract <- function(variants, bfile, plink_bin) {
  variants <- unique(as.character(variants))
  if (length(variants) == 0) return(NULL)

  fn <- tempfile()
  on.exit(unlink(paste0(fn, "*")), add = TRUE)

  tryCatch({
    data.table::fwrite(list(variants), fn, col.names = FALSE)

    if (!file.exists(paste0(bfile, ".bed"))) {
      message("[LD] PLINK bfile not found: ", bfile)
      return(NULL)
    }

    shell <- ifelse(Sys.info()["sysname"] == "Darwin", "cmd", "sh")
    cmd <- paste(
      shQuote(plink_bin, type = shell),
      "--bfile", shQuote(bfile, type = shell),
      "--extract", shQuote(fn, type = shell),
      "--r square",
      "--make-just-bim",
      "--out", shQuote(fn, type = shell)
    )

    sys_out <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)

    ld_file <- paste0(fn, ".ld")
    bim_file <- paste0(fn, ".bim")

    if (!file.exists(ld_file)) {
      message("[LD] LD file not created. Check PLINK output.")
      return(NULL)
    }

    ld_mat <- as.matrix(data.table::fread(ld_file))
    if (!file.exists(bim_file)) {
      message("[LD] BIM file not found")
      return(NULL)
    }
    bim <- data.table::fread(bim_file, col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"))

    colnames(ld_mat) <- bim$SNP
    rownames(ld_mat) <- bim$SNP

    ld_result <- as.data.frame(as.table(ld_mat))
    colnames(ld_result) <- c("SNP_A", "SNP_B", "R")
    ld_result <- ld_result[ld_result$SNP_A != ld_result$SNP_B, ]

    n_snps_in_bim <- nrow(bim)
    message(glue("[LD] Extracted LD for {length(unique(c(ld_result$SNP_A, ld_result$SNP_B)))}/{length(variants)} SNPs (PLINK has {n_snps_in_bim} SNPs in region)"))

    res_inv <- data.frame(SNP_A = ld_result$SNP_B, SNP_B = ld_result$SNP_A, R = ld_result$R)

    return(dplyr::bind_rows(ld_result, res_inv) %>% dplyr::distinct())
  }, error = function(e) {
    message(glue("[LD] Extraction failed: {e$message}"))
    return(NULL)
  })
}

LD_plot <- function(df, ld_df = NULL, lead_snps = NULL,
                    bfile = NULL, plink_bin = "plink",
                    plot_title = NULL, plot_subtitle = NULL,
                    region_recomb = NULL, xlim = NULL,
                    show_lead_line = FALSE) {

   df <- safe_df_convert(df)
  if (is.null(df)) {
    return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No Input Data")))
  }

  col_snp <- resolve_col(df, "rsid", c("SNPID.gwas", "snp", "SNP", "SNPID", "marker"))
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
  rsid_sample <- head(plot_df$rsid[!is.na(plot_df$rsid)], 3)
  is_rsid_format <- all(grepl("^rs", rsid_sample))
  message(glue("[LD_plot] SNP column: {col_snp}, Sample: {paste(rsid_sample, collapse=', ')}, Is rsID format: {is_rsid_format}"))

  plot_df <- plot_df[
    !is.na(plot_df$rsid) & plot_df$rsid != "" &
    !is.na(plot_df$position) &
    !is.na(plot_df$p_value) & plot_df$p_value > 0,
  ]

  if (nrow(plot_df) == 0) {
    return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No Valid SNPs")))
  }

  if (!is.null(xlim) && length(xlim) == 2) {
    plot_df <- plot_df[plot_df$position >= xlim[1] & plot_df$position <= xlim[2], ]
    if (nrow(plot_df) == 0) {
      return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No SNPs in Window")))
    }
  }

  message(glue("[LD_plot] Plotting {nrow(plot_df)} SNPs"))
  if (!is.null(lead_snps) && length(lead_snps) > 0) {
    lead_snp <- lead_snps[1]
    if (!lead_snp %in% plot_df$rsid) {
      lead_snp <- plot_df$rsid[which.min(plot_df$p_value)]
    }
  } else {
    lead_snp <- plot_df$rsid[which.min(plot_df$p_value)]
  }
  lead_pos <- plot_df$position[plot_df$rsid == lead_snp][1]

   message(glue("[LD_plot] Lead SNP: {lead_snp}"))
   message(glue("[LD_plot] Lead SNP in data: {lead_snp %in% plot_df$rsid}"))
   if (!(lead_snp %in% plot_df$rsid)) {
     close_matches <- grep(lead_snp, plot_df$rsid, value = TRUE)
     message(glue("[LD] Lead SNP '{lead_snp}' not found. Close matches: {paste(head(close_matches, 5), collapse=', ')}"))
     trimmed_match <- grep(gsub("^\\s+|\\s+$", "", lead_snp), plot_df$rsid, value = TRUE)
     if (length(trimmed_match) > 0) {
       message(glue("[LD] Found after trimming whitespace"))
     }
   }
   sample_snps <- head(unique(plot_df$rsid), 5)
   message(glue("[LD_plot] Sample SNPs in data: {paste(sample_snps, collapse=', ')}"))

   plot_df$r2 <- NA_real_
   if (!is.null(bfile) && file.exists(paste0(bfile, ".bed"))) {
     snps_for_ld <- unique(c(lead_snp, plot_df$rsid))
     message(glue("[LD] Extracting LD for {length(snps_for_ld)} SNPs including lead SNP..."))
     
     ld_result <- ld_extract(snps_for_ld, bfile, plink_bin)

     if (!is.null(ld_result)) {
       lead_in_ld <- any(ld_result$SNP_A == lead_snp | ld_result$SNP_B == lead_snp)
       message(glue("[LD] Lead SNP in LD results: {lead_in_ld}"))
       
       if (lead_in_ld) {
         ld_to_lead <- ld_result[ld_result$SNP_A == lead_snp, ]
         message(glue("[LD_plot] LD pairs with lead SNP: {nrow(ld_to_lead)}"))
        } else {
          message(glue("[LD] WARNING: Lead SNP {lead_snp} not in PLINK reference!"))
          snps_in_ld <- unique(c(ld_result$SNP_A, ld_result$SNP_B))
          plot_df$in_ld <- plot_df$rsid %in% snps_in_ld
          
          if (any(plot_df$in_ld)) {
            snps_in_ld_df <- plot_df[plot_df$in_ld, ]
            if (nrow(snps_in_ld_df) > 0) {
              alt_lead <- snps_in_ld_df$rsid[which.min(snps_in_ld_df$p_value)]
              message(glue("[LD] Using {alt_lead} as alternative lead SNP for LD visualization"))
              ld_to_lead <- ld_result[ld_result$SNP_A == alt_lead, ]
              if (nrow(ld_to_lead) > 0) {
                ld_to_lead$r2_calc <- abs(ld_to_lead$R)^2
                plot_df <- merge(
                  plot_df,
                  ld_to_lead[, c("SNP_B", "r2_calc")],
                  by.x = "rsid",
                  by.y = "SNP_B",
                  all.x = TRUE
                )
                
                plot_df$r2 <- ifelse(!is.na(plot_df$r2_calc), plot_df$r2_calc, plot_df$r2)
                plot_df$r2_calc <- NULL
                plot_df$r2[plot_df$rsid == alt_lead] <- 1.0
                lead_snp <- alt_lead
              }
            }
          }
          ld_to_lead <- data.frame()
        }

       if (nrow(ld_to_lead) > 0) {
         ld_to_lead$r2_calc <- abs(ld_to_lead$R)^2

         plot_df <- merge(
           plot_df,
           ld_to_lead[, c("SNP_B", "r2_calc")],
           by.x = "rsid",
           by.y = "SNP_B",
           all.x = TRUE
         )

         plot_df$r2 <- ifelse(!is.na(plot_df$r2_calc), plot_df$r2_calc, plot_df$r2)
         plot_df$r2_calc <- NULL
       }
     } else {
       message("[LD_plot] LD extraction returned NULL - check PLINK bfile and SNP IDs")
     }
   } else {
     message(glue("[LD_plot] PLINK bfile not found: {bfile}"))
   }
   plot_df$r2[plot_df$rsid == lead_snp] <- 1.0
   plot_df$r2[is.na(plot_df$r2)] <- 0.05

   high_ld_count <- sum(plot_df$r2 > 0.2 & plot_df$rsid != lead_snp)
   message(glue("[LD_plot] Found {high_ld_count} proxy SNPs with r² > 0.2"))

   plot_df$is_lead <- plot_df$rsid == lead_snp
   plot_df$color_code <- cut(
    plot_df$r2,
    breaks = c(-0.01, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c("#313695", "#4575B4", "#74ADD1", "#FDB863", "#D73027"),
    include.lowest = TRUE
  )
  plot_df$color_code <- as.character(plot_df$color_code)
  plot_df$color_code[plot_df$is_lead] <- "#7F3C8D"

  maxlogP <- max(-log10(plot_df$p_value), na.rm = TRUE)
  if (is.infinite(maxlogP) || maxlogP < 1) maxlogP <- 10
  p <- ggplot()
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
          scale_factor <- (maxlogP * 0.85) / recomb_max
          region_recomb$scaled_rate <- region_recomb$rate_mapped * scale_factor

          p <- p +
            geom_ribbon(data = region_recomb,
                        aes(x = pos_mapped/1e6, ymin = 0, ymax = scaled_rate),
                        fill = "#DEEBF7", alpha = 0.5, inherit.aes = FALSE) +
            geom_line(data = region_recomb,
                      aes(x = pos_mapped/1e6, y = scaled_rate),
                      color = "#3182BD", linewidth = 0.6, alpha = 0.8, inherit.aes = FALSE)
        }
      }
    }
  }
  if (show_lead_line && !is.na(lead_pos)) {
    p <- p + 
      geom_vline(xintercept = lead_pos/1e6, 
                 linetype = "dotted", 
                 color = "#7F3C8D", 
                 linewidth = 0.8, 
                 alpha = 0.6)
  }
  p <- p +
    geom_point(data = plot_df,
               aes(x = position/1e6, y = -log10(p_value), 
                   fill = color_code, size = is_lead, shape = is_lead),
               alpha = 0.9, stroke = 0.4, color = "white") +
    geom_hline(yintercept = -log10(5e-8), 
               linetype = "dashed", color = "grey30", linewidth = 0.6) +
    scale_fill_identity(
      name = expression(italic(r)^2 ~ "with lead SNP"),
      labels = c("Lead SNP", "0.8 – 1.0", "0.6 – 0.8", "0.4 – 0.6", "0.2 – 0.4", "< 0.2"),
      breaks = c("#7F3C8D", "#D73027", "#FDB863", "#74ADD1", "#4575B4", "#313695"),
      guide = guide_legend(
        override.aes = list(shape = 21, size = 4, alpha = 1, stroke = 0.4, color = "white"),
        title.position = "top",
        title.hjust = 0.5,
        nrow = 6,
        byrow = TRUE
      )
    ) +
    scale_size_manual(values = c(2.5, 6), guide = "none") +
    scale_shape_manual(values = c(21, 23), guide = "none") +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = NULL,
      y = expression(-log[10] ~ italic(P))
    ) +
    theme_classic(base_size = 11, base_family = "sans") +
    theme(
      legend.position = "right",
      legend.justification = c(0, 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.7, "lines"),
      legend.key.height = unit(0.8, "lines"),
      legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
      legend.margin = margin(t = 3, r = 3, b = 3, l = 3),
      legend.box.spacing = unit(0, "pt"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.line.x = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.12, "cm"),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 6)),
      plot.title = element_text(size = 12, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0, margin = margin(b = 4)),
      plot.margin = margin(t = 8, r = 5, b = 2, l = 8, unit = "pt"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  if (!is.null(xlim) && length(xlim) == 2) {
    p <- p + scale_x_continuous(limits = xlim/1e6, expand = expansion(mult = c(0.01, 0.01)))
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)))
  }
  if (!is.null(region_recomb) && exists("scale_factor")) {
    p <- p + scale_y_continuous(
      name = expression(-log[10] ~ italic(P)),
      expand = expansion(mult = c(0, 0.05)),
      sec.axis = sec_axis(
        ~ . / scale_factor,
        name = "Recombination rate (cM/Mb)"
      )
    ) +
    theme(
      axis.title.y.right = element_text(size = 10, face = "bold", margin = margin(l = 6)),
      axis.text.y.right = element_text(size = 9, color = "black")
    )
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  }

  return(p)
}

genetrack <- function(chrom_str, xlim_bp, gtf_path = NULL, show_strand_legend = FALSE) {
    chrom_clean <- gsub("^chr", "", as.character(chrom_str))
    chrom_str <- paste0("chr", chrom_clean)
    
    if (is.null(xlim_bp) || length(xlim_bp) != 2) {
      xlim_bp <- c(0, 1e6)
    }

    BPStart <- xlim_bp[1]
    BPStop <- xlim_bp[2]

    gr <- GRanges(seqnames = chrom_str, ranges = IRanges(start = BPStart, end = BPStop))
    tx_df <- NULL
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
                
                if ("strand" %in% names(df_gtf)) {
                    df_gtf$gene_strand <- as.character(df_gtf$strand)
                } else {
                    df_gtf$gene_strand <- "+"
                }
                
                df_gtf$gene_strand[is.na(df_gtf$gene_strand) | df_gtf$gene_strand == "*"] <- "+"

                if ("gene_type" %in% names(df_gtf)) {
                    df_gtf <- df_gtf[df_gtf$gene_type == "protein_coding", ]
                }

                if (nrow(df_gtf) > 0) {
                    tx_summary <- df_gtf %>%
                        dplyr::group_by(symbol) %>%
                        dplyr::summarize(
                            start = min(start, na.rm = TRUE),
                            end = max(end, na.rm = TRUE),
                            strand = dplyr::first(gene_strand),
                            .groups = 'drop'
                        )
                    tx_df <- as.data.frame(tx_summary)
                }
            }
        }, error = function(e) {
            message(glue("[Gene Track] GTF error: {e$message}"))
        })
    }
    if (is.null(tx_df) || nrow(tx_df) == 0) {
        tryCatch({
            txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
            subset_tx <- GenomicFeatures::transcriptsByOverlaps(txdb, gr)

            if (length(subset_tx) > 0) {
                gene_ids <- subset_tx$gene_id
                gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
                
                symbols <- NULL
                if (length(gene_ids) > 0) {
                    symbols <- suppressMessages(
                        AnnotationDbi::mapIds(
                            org.Hs.eg.db::org.Hs.eg.db,
                            keys = gene_ids,
                            column = "SYMBOL",
                            keytype = "ENTREZID",
                            multiVals = "first"
                        )
                    )
                }

                tx_df <- data.frame(
                    gene_id = subset_tx$gene_id,
                    tx_name = subset_tx$tx_name,
                    start = start(subset_tx),
                    end = end(subset_tx),
                    strand = as.character(strand(subset_tx)),
                    stringsAsFactors = FALSE
                )
                
                if (!is.null(symbols) && length(symbols) > 0) {
                    symbol_df <- data.frame(
                        gene_id = names(symbols),
                        symbol = as.character(symbols),
                        stringsAsFactors = FALSE
                    )
                    tx_df <- merge(tx_df, symbol_df, by = "gene_id", all.x = TRUE)
                    tx_df$symbol <- ifelse(is.na(tx_df$symbol), tx_df$tx_name, tx_df$symbol)
                } else {
                    tx_df$symbol <- tx_df$tx_name
                }

                tx_df <- tx_df[!grepl("^(LOC|LINC|MIR|SNOR|RN7|RNU)", tx_df$symbol), ]

                if (nrow(tx_df) > 0) {
                    tx_df <- tx_df %>%
                        dplyr::group_by(symbol) %>%
                        dplyr::summarize(
                          start = min(start, na.rm = TRUE), 
                          end = max(end, na.rm = TRUE),
                          strand = dplyr::first(strand),
                          .groups = 'drop'
                        ) %>%
                        as.data.frame()
                }
            }
        }, error = function(e) {
            message(glue("[Gene Track] TxDb error: {e$message}"))
        })
    }

    if (is.null(tx_df) || nrow(tx_df) == 0) {
        return(
          ggplot() + 
          theme_void() +
          geom_text(aes(x = mean(xlim_bp)/1e6, y = 0.5, label = "No genes in region"), 
                    size = 3, color = "grey50", fontface = "italic") +
          scale_x_continuous(
            limits = xlim_bp/1e6, 
            expand = c(0.01, 0.01),
            name = paste0("Position on ", chrom_str, " (Mb)")
          ) +
          theme(
            plot.margin = margin(t = 2, r = 5, b = 8, l = 8, unit = "pt"),
            axis.text.x = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 6)),
            axis.ticks.x = element_line(color = "black", linewidth = 0.5),
            axis.ticks.length.x = unit(0.12, "cm"),
            axis.line.x = element_line(color = "black", linewidth = 0.6)
          )
        )
    }

    message(glue("[Gene Track] Drawing {nrow(tx_df)} genes"))

    tx_df <- tx_df[order(tx_df$start), ]
    
    tx_df$strand[is.na(tx_df$strand) | tx_df$strand == "*" | tx_df$strand == ""] <- "+"
    tx_df$strand <- ifelse(tx_df$strand %in% c("+", "-"), tx_df$strand, "+")
    tx_df$y <- 1
    if (nrow(tx_df) > 1) {
        for (i in 2:nrow(tx_df)) {
            overlap_margin <- (BPStop - BPStart) * 0.03
            if (tx_df$start[i] < (tx_df$end[i-1] + overlap_margin)) {
                tx_df$y[i] <- ifelse(tx_df$y[i-1] == 1, 2, 1)
            }
        }
    }
    
    tx_df$color <- ifelse(tx_df$strand == "+", "#D73027", "#4575B4")
    p <- ggplot(tx_df)
    p <- p + geom_segment(
      aes(x = start/1e6, xend = end/1e6, y = y, yend = y, color = strand),
      linewidth = 2,
      lineend = "round"
    )
    arrow_length <- (BPStop - BPStart) * 0.015
    
    fwd_genes <- tx_df[tx_df$strand == "+", ]
    if (nrow(fwd_genes) > 0) {
      p <- p + geom_segment(
        data = fwd_genes,
        aes(x = (end - arrow_length)/1e6, xend = end/1e6, y = y, yend = y),
        arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
        color = "#D73027",
        linewidth = 1.1
      )
    }
    
    rev_genes <- tx_df[tx_df$strand == "-", ]
    if (nrow(rev_genes) > 0) {
      p <- p + geom_segment(
        data = rev_genes,
        aes(x = (start + arrow_length)/1e6, xend = start/1e6, y = y, yend = y),
        arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
        color = "#4575B4",
        linewidth = 1.1
      )
    }
    p <- p + geom_text(
      aes(x = (start + end)/2e6, y = y, label = symbol, color = strand),
      size = 3, 
      fontface = "bold.italic",
      vjust = -1.1,
      show.legend = FALSE
    )
    p <- p +
      scale_color_manual(
        values = c("+" = "#D73027", "-" = "#4575B4"),
        labels = c("+" = "Forward", "-" = "Reverse"),
        name = "Strand"
      ) +
      scale_x_continuous(
        limits = xlim_bp/1e6, 
        expand = expansion(mult = c(0.01, 0.01)),
        name = paste0("Position on ", chrom_str, " (Mb)")
      ) +
      scale_y_continuous(limits = c(0.3, 2.7), expand = c(0, 0)) +
      theme_void() +
      theme(
        plot.margin = margin(t = 2, r = 5, b = 8, l = 8, unit = "pt"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 6)),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length.x = unit(0.12, "cm"),
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        legend.position = if(show_strand_legend) "bottom" else "none",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8)
      )

    return(p)
}

plot_qtl_association <- function(qtl_all_chrom, qtl_all_pvalue, leadSNP_DF,
                                 ld_df = NULL, gtf_path = NULL, region_recomb = NULL,
                                 show_lead_line = FALSE) {

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
            plot_title = paste(lead_snp_val, "–", gene_sym, "QTL"),
            plot_subtitle = "eQTL Association (hg38)",
            region_recomb = region_recomb,
            xlim = c(min_pos, max_pos),
            show_lead_line = show_lead_line
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
    ggarrange(p_assoc, p_track, 
              ncol = 1, 
              heights = c(3.5, 1.5), 
              align = "v") 
}

plot_gwas_association <- function(colocInputFile, qtl_all_chrom, leadSNP_DF,
                                  ld_df = NULL, gtf_path = NULL, region_recomb = NULL,
                                  show_lead_line = FALSE) {

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
            plot_title = paste(lead_snp_val, "–", trait_name, "GWAS"),
            plot_subtitle = "GWAS Association (hg38)",
            region_recomb = region_recomb,
            xlim = c(min_pos, max_pos),
            show_lead_line = show_lead_line
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
    ggarrange(p_assoc, p_track, 
              ncol = 1, 
              heights = c(3.5, 1.5),
              align = "v")
}
