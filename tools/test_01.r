# ------------------------------------------------------------------------------
# src/utils_plot.R - Ultra Robust Version
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
# Helper Functions (Enhanced)
# ==============================================================================

#' Resolve column name with fallback
resolve_col <- function(df, preferred, alternatives = NULL) {
  if (preferred %in% names(df)) return(preferred)
  if (!is.null(alternatives)) {
    for (alt in alternatives) {
      if (alt %in% names(df)) return(alt)
    }
  }
  # 如果都找不到，返回 NULL 而不是第一列
  warning(glue("Column '{preferred}' not found. Alternatives: {paste(alternatives, collapse=', ')}"))
  return(NULL)
}

#' Safe data.frame conversion with validation
safe_df_convert <- function(df) {
    if (is.null(df)) {
        message("[DEBUG] Input is NULL")
        return(NULL)
    }
    
    # 转换为纯 data.frame
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    
    if (nrow(df) == 0) {
        message("[DEBUG] Input has 0 rows")
        return(NULL)
    }
    
    # 打印列名用于调试
    message(glue("[DEBUG] Available columns: {paste(names(df), collapse=', ')}"))
    
    return(df)
}

#' Extract numeric vector safely
safe_numeric <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0) return(default)
    result <- suppressWarnings(as.numeric(as.character(x)))
    result[is.na(result)] <- default
    return(result)
}

#' Extract character vector safely
safe_character <- function(x, default = NA_character_) {
    if (is.null(x) || length(x) == 0) return(default)
    result <- as.character(x)
    result[is.na(result) | result == ""] <- default
    return(result)
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
# LD Plot (Simplified & Robust)
# ==============================================================================

LD_plot <- function(df, ld_df = NULL, lead_snps = NULL,
                    bfile = NULL, plink_bin = "plink",
                    plot_title = NULL, plot_subtitle = NULL,
                    region_recomb = NULL) {
  
  # 数据验证
  df <- safe_df_convert(df)
  if (is.null(df)) {
    return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No Input Data")))
  }
  
  # 查找必需列
  col_snp <- resolve_col(df, "rsid", c("snp", "SNP", "SNPID", "marker"))
  col_chr <- resolve_col(df, "chromosome", c("CHR", "chr", "CHR.qtl", "CHR.gwas"))
  col_pos <- resolve_col(df, "position", c("POS", "pos", "POS.qtl", "POS.gwas", "BP"))
  col_p   <- resolve_col(df, "p_value", c("P", "pval", "P.qtl", "P.gwas"))
  
  # 验证所有必需列存在
  if (is.null(col_snp) || is.null(col_chr) || is.null(col_pos) || is.null(col_p)) {
    missing <- c(
      if(is.null(col_snp)) "SNP ID",
      if(is.null(col_chr)) "Chromosome",
      if(is.null(col_pos)) "Position",
      if(is.null(col_p)) "P-value"
    )
    msg <- paste("Missing columns:", paste(missing, collapse=", "))
    message(glue("[LD_plot ERROR] {msg}"))
    return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label=msg)))
  }
  
  # 构建干净的数据框
  plot_df <- data.frame(
    rsid = safe_character(df[[col_snp]]),
    chromosome = safe_character(df[[col_chr]]),
    position = safe_numeric(df[[col_pos]]),
    p_value = safe_numeric(df[[col_p]]),
    stringsAsFactors = FALSE
  )
  
  # 移除无效行
  plot_df <- plot_df[
    !is.na(plot_df$rsid) & plot_df$rsid != "" &
    !is.na(plot_df$position) & 
    !is.na(plot_df$p_value) & plot_df$p_value > 0,
  ]
  
  if (nrow(plot_df) == 0) {
    return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No Valid SNPs")))
  }
  
  message(glue("[LD_plot] Plotting {nrow(plot_df)} SNPs"))
  
  # 确定 lead SNP
  if (!is.null(lead_snps) && length(lead_snps) > 0) {
    lead_snp <- lead_snps[1]
    if (!lead_snp %in% plot_df$rsid) {
      message(glue("[LD_plot] Lead SNP '{lead_snp}' not in data, using most significant"))
      lead_snp <- plot_df$rsid[which.min(plot_df$p_value)]
    }
  } else {
    lead_snp <- plot_df$rsid[which.min(plot_df$p_value)]
  }
  
  message(glue("[LD_plot] Lead SNP: {lead_snp}"))
  
  # 计算 LD（如果可能）
  plot_df$r2 <- 0.1  # 默认值
  
  if (!is.null(bfile) && file.exists(paste0(bfile, ".bed"))) {
    ld_result <- ld_extract(unique(c(lead_snp, plot_df$rsid)), bfile, plink_bin)
    
    if (!is.null(ld_result)) {
      # 只保留与 lead SNP 的 LD
      ld_to_lead <- ld_result[ld_result$SNP_A == lead_snp, ]
      if (nrow(ld_to_lead) > 0) {
        ld_to_lead$r2 <- abs(ld_to_lead$R)^2
        plot_df <- merge(
          plot_df, 
          ld_to_lead[, c("SNP_B", "r2")], 
          by.x = "rsid", 
          by.y = "SNP_B", 
          all.x = TRUE,
          suffixes = c("", "_ld")
        )
        plot_df$r2[is.na(plot_df$r2)] <- 0.1
      }
    }
  }
  
  # 标记 lead SNP
  plot_df$is_lead <- plot_df$rsid == lead_snp
  
  # 颜色编码
  plot_df$color_code <- cut(
    plot_df$r2, 
    breaks = c(-0.1, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c("blue4", "blue", "darkgreen", "orange", "red"),
    include.lowest = TRUE
  )
  plot_df$color_code <- as.character(plot_df$color_code)
  plot_df$color_code[plot_df$is_lead] <- "purple"
  
  # 绘图
  maxlogP <- max(-log10(plot_df$p_value), na.rm = TRUE)
  if (is.infinite(maxlogP) || maxlogP < 1) maxlogP <- 10
  
  p <- ggplot(plot_df, aes(x = position/1e6, y = -log10(p_value))) +
    geom_point(aes(fill = color_code, size = is_lead, shape = is_lead), alpha = 0.8) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey") +
    scale_fill_identity(
      name = expression(LD ~ (r^2)),
      labels = c("Lead SNP", "0.8-1.0", "0.6-0.8", "0.4-0.6", "0.2-0.4", "0-0.2"),
      breaks = c("purple", "red", "orange", "darkgreen", "blue", "blue4"),
      guide = guide_legend(override.aes = list(shape = 21, size = 5, alpha = 1))
    ) +
    scale_size_manual(values = c(3, 6), guide = "none") +
    scale_shape_manual(values = c(21, 23), guide = "none") +
    scale_y_continuous(name = expression(-log[10]("p-value"))) +
    labs(title = plot_title, subtitle = plot_subtitle, x = "Position (Mb)") +
    theme_light(base_size = 12) +
    theme(legend.position = "right")
  
  # 添加重组率（如果有）
  if (!is.null(region_recomb) && is.data.frame(region_recomb) && nrow(region_recomb) > 0) {
    if ("Rate(cM/Mb)" %in% names(region_recomb) && "Position(bp)" %in% names(region_recomb)) {
      recomb_max <- max(region_recomb$`Rate(cM/Mb)`, na.rm = TRUE)
      if (!is.infinite(recomb_max) && recomb_max > 0) {
        scale_factor <- maxlogP / recomb_max
        region_recomb$scaled_rate <- region_recomb$`Rate(cM/Mb)` * scale_factor
        
        p <- p +
          geom_line(
            data = region_recomb,
            aes(x = `Position(bp)`/1e6, y = scaled_rate),
            color = "blue", linewidth = 0.8, alpha = 0.6, inherit.aes = FALSE
          ) +
          scale_y_continuous(
            name = expression(-log[10]("p-value")),
            sec.axis = sec_axis(~ . / scale_factor, name = "Recomb Rate (cM/Mb)")
          )
      }
    }
  }
  
  return(p)
}

# ==============================================================================
# Gene Track (Simplified)
# ==============================================================================

genetrack <- function(chrom_str, BPStart, BPStop, gtf_path = NULL) {
    chrom_clean <- gsub("^chr", "", as.character(chrom_str))
    chrom_str <- paste0("chr", chrom_clean)
    
    gr <- GRanges(seqnames = chrom_str, ranges = IRanges(start = BPStart, end = BPStop))
    
    tx_df <- NULL
    
    # 尝试从 GTF 加载
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
                
                # 每个基因只保留一个转录本
                tx_summary <- df_gtf %>%
                    group_by(symbol) %>%
                    summarize(
                        start = min(start),
                        end = max(end),
                        .groups = 'drop'
                    )
                
                tx_df <- as.data.frame(tx_summary)
            }
        }, error = function(e) {
            message(glue("[Gene Track] GTF error: {e$message}"))
        })
    }
    
    # 回退到 TxDb
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
                
                # 每个基因一行
                tx_df <- tx_df %>%
                    group_by(symbol) %>%
                    summarize(start = min(start), end = max(end), .groups = 'drop') %>%
                    as.data.frame()
            }
        }, error = function(e) {
            message(glue("[Gene Track] TxDb error: {e$message}"))
        })
    }
    
    # 如果还是没有数据，返回空图
    if (is.null(tx_df) || nrow(tx_df) == 0) {
        return(ggplot() + theme_void() + 
               geom_text(aes(x = (BPStart+BPStop)/2, y = 0.5, label = "No genes")))
    }
    
    # 绘制基因
    tx_df <- tx_df[order(tx_df$start), ]
    tx_df$y <- seq_len(nrow(tx_df))
    
    ggplot(tx_df) +
        geom_segment(aes(x = start, xend = end, y = y, yend = y), 
                     color = "darkblue", linewidth = 1) +
        geom_text(aes(x = (start + end)/2, y = y, label = symbol),
                  size = 3, fontface = "italic", vjust = -0.5) +
        scale_x_continuous(limits = c(BPStart, BPStop), expand = c(0, 0)) +
        theme_void() +
        theme(plot.margin = margin(t = 0, unit = "cm"))
}

# ==============================================================================
# Composite Functions
# ==============================================================================

plot_qtl_association <- function(qtl_all_chrom, qtl_all_pvalue, leadSNP_DF,
                                 ld_df = NULL, gtf_path = NULL, region_recomb = NULL) {
    
    message("[plot_qtl_association] Starting...")
    
    leadSNP_DF <- safe_df_convert(leadSNP_DF)
    if (is.null(leadSNP_DF)) {
        return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No data")))
    }
    
    # 获取全局变量
    lead_snp_val <- if(exists("lead_SNP", envir=.GlobalEnv)) get("lead_SNP", envir=.GlobalEnv) else NULL
    gene_sym <- if(exists("geneSymbol", envir=.GlobalEnv)) get("geneSymbol", envir=.GlobalEnv) else "Gene"
    plink_bfile_val <- if(exists("plink_bfile", envir=.GlobalEnv)) get("plink_bfile", envir=.GlobalEnv) else NULL
    
    # 确定绘图区域（基于所有 SNP）
    pos_col <- resolve_col(leadSNP_DF, "POS.qtl", c("POS", "pos", "BP", "POS.gwas"))
    chr_col <- resolve_col(leadSNP_DF, qtl_all_chrom, c("CHR", "CHR.qtl", "chr"))
    
    if (is.null(pos_col) || is.null(chr_col)) {
        return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="Missing position columns")))
    }
    
    positions <- safe_numeric(leadSNP_DF[[pos_col]])
    positions <- positions[!is.na(positions) & positions > 0]
    
    if (length(positions) == 0) {
        return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No valid positions")))
    }
    
    min_pos <- min(positions)
    max_pos <- max(positions)
    chrom_num <- safe_character(leadSNP_DF[[chr_col]])[1]
    
    message(glue("[plot_qtl_association] Region: chr{chrom_num}:{min_pos}-{max_pos}"))
    
    # 生成关联图
    p_assoc <- tryCatch({
        LD_plot(
            df = leadSNP_DF,
            ld_df = ld_df,
            lead_snps = lead_snp_val,
            bfile = plink_bfile_val,
            plink_bin = "plink",
            plot_title = paste(lead_snp_val, gene_sym, "QTL"),
            plot_subtitle = "QTL Association (hg38)",
            region_recomb = region_recomb
        )
    }, error = function(e) {
        message(glue("[plot_qtl_association] LD_plot error: {e$message}"))
        ggplot() + theme_void() + geom_text(aes(x=0, y=0, label=paste("LD plot failed:", e$message)))
    })
    
    # 生成基因轨迹
    p_track <- tryCatch({
        genetrack(chrom_num, min_pos, max_pos, gtf_path = gtf_path)
    }, error = function(e) {
        message(glue("[plot_qtl_association] Gene track error: {e$message}"))
        ggplot() + theme_void()
    })
    
    # 组合
    ggarrange(p_assoc, p_track, ncol = 1, heights = c(3, 1), align = "v")
}

plot_gwas_association <- function(colocInputFile, qtl_all_chrom, leadSNP_DF,
                                  ld_df = NULL, gtf_path = NULL, region_recomb = NULL) {
    
    message("[plot_gwas_association] Starting...")
    
    leadSNP_DF <- safe_df_convert(leadSNP_DF)
    if (is.null(leadSNP_DF)) {
        return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No data")))
    }
    
    # 获取全局变量
    lead_snp_val <- if(exists("lead_SNP", envir=.GlobalEnv)) get("lead_SNP", envir=.GlobalEnv) else NULL
    trait_name <- if(exists("trait", envir=.GlobalEnv)) get("trait", envir=.GlobalEnv) else "GWAS"
    plink_bfile_val <- if(exists("plink_bfile", envir=.GlobalEnv)) get("plink_bfile", envir=.GlobalEnv) else NULL
    
    # 确定绘图区域（使用 GWAS 数据的列）
    pos_col <- resolve_col(leadSNP_DF, "POS.gwas", c("POS", "pos", "BP", "POS.qtl"))
    chr_col <- resolve_col(leadSNP_DF, "CHR.gwas", c("CHR", "chr", "CHR.qtl"))
    
    if (is.null(pos_col) || is.null(chr_col)) {
        return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="Missing position columns")))
    }
    
    positions <- safe_numeric(leadSNP_DF[[pos_col]])
    positions <- positions[!is.na(positions) & positions > 0]
    
    if (length(positions) == 0) {
        return(ggplot() + theme_void() + geom_text(aes(x=0, y=0, label="No valid positions")))
    }
    
    min_pos <- min(positions)
    max_pos <- max(positions)
    chrom_num <- safe_character(leadSNP_DF[[chr_col]])[1]
    
    message(glue("[plot_gwas_association] Region: chr{chrom_num}:{min_pos}-{max_pos}"))
    
    # 生成关联图（使用 GWAS p-value）
    p_assoc <- tryCatch({
        LD_plot(
            df = leadSNP_DF,
            ld_df = ld_df,
            lead_snps = lead_snp_val,
            bfile = plink_bfile_val,
            plink_bin = "plink",
            plot_title = paste(lead_snp_val, trait_name, "GWAS"),
            plot_subtitle = "GWAS Association (hg38)",
            region_recomb = region_recomb
        )
    }, error = function(e) {
        message(glue("[plot_gwas_association] LD_plot error: {e$message}"))
        ggplot() + theme_void() + geom_text(aes(x=0, y=0, label=paste("LD plot failed:", e$message)))
    })
    
    # 生成基因轨迹
    p_track <- tryCatch({
        genetrack(chrom_num, min_pos, max_pos, gtf_path = gtf_path)
    }, error = function(e) {
        message(glue("[plot_gwas_association] Gene track error: {e$message}"))
        ggplot() + theme_void()
    })
    
    # 组合
    ggarrange(p_assoc, p_track, ncol = 1, heights = c(3, 1), align = "v")
}
