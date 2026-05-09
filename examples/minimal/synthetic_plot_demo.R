#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

make_synthetic_locus <- function() {
  set.seed(20260327)

  n_bg <- 80
  bg_pos <- sort(sample(seq(43120000, 43580000, by = 500), n_bg))
  bg_p <- runif(n_bg, min = 5e-3, max = 0.08)

  signal_pos <- c(43562000, 43564000, 43566000, 43568000, 43570000, 43572000)
  signal_p <- c(2e-6, 8e-7, 2e-8, 7e-8, 3e-6, 9e-6)
  signal_rs <- c("rs2819340", "rs3791138", "rs10890255", "rs3791137", "rs2004899", "rs2819341")

  df <- data.frame(
    rsid = c(paste0("rsBG", seq_len(n_bg)), signal_rs),
    CHR.qtl = "chr1",
    POS.qtl = c(bg_pos, signal_pos),
    P.qtl = c(bg_p, signal_p),
    P.gwas = c(runif(n_bg, min = 1e-3, max = 0.2), c(3e-5, 6e-6, 1e-9, 5e-6, 9e-5, 7e-5)),
    stringsAsFactors = FALSE
  )

  df$category <- "Other variants"
  df$category[df$rsid %in% c("rs2819340", "rs3791138", "rs2004899")] <- "95% CS"
  df$category[df$rsid == "rs10890255"] <- "Lead SNP"
  df$category <- factor(df$category, levels = c("Other variants", "95% CS", "Lead SNP"))
  df$mb <- df$POS.qtl / 1e6
  df$logp <- -log10(df$P.qtl)
  df
}

synthetic_df <- make_synthetic_locus()

palette <- c(
  "Other variants" = "#8A8F98",
  "95% CS" = "#0072B2",
  "Lead SNP" = "#D55E00"
)

plot_obj <- ggplot(synthetic_df, aes(x = mb, y = logp)) +
  geom_hline(yintercept = -log10(5e-6), linewidth = 0.25, linetype = "22", color = "#B8B8B8") +
  annotate("segment", x = 43.46, xend = 43.51, y = -0.38, yend = -0.38,
           linewidth = 1.6, lineend = "round", color = "#222222") +
  annotate("text", x = 43.512, y = -0.38, label = "CCDC30", hjust = 0, size = 1.75, color = "#222222") +
  geom_point(
    aes(fill = category, size = category),
    shape = 21,
    stroke = 0.25,
    color = "white",
    alpha = 0.96
  ) +
  scale_fill_manual(values = palette, name = NULL) +
  scale_size_manual(values = c("Other variants" = 1.5, "95% CS" = 2.2, "Lead SNP" = 2.8), guide = "none") +
  scale_x_continuous(
    breaks = c(43.2, 43.4, 43.6),
    labels = c("43.2", "43.4", "43.6"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(-0.58, 8.2),
    breaks = c(0, 2, 4, 6, 8),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "CCDC30 locus",
    x = "Position on chr1 (Mb)",
    y = expression(-log[10](QTL~italic(P)))
  ) +
  guides(fill = guide_legend(override.aes = list(size = c(2.2, 2.6, 3.0), alpha = 1))) +
  theme_classic(base_size = 6.8) +
  theme(
    plot.title = element_text(size = 7.8, face = "bold", margin = margin(b = 2)),
    axis.title = element_text(size = 6.8),
    axis.text = element_text(size = 6.0, color = "#333333"),
    axis.line = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks = element_line(linewidth = 0.25, color = "#333333"),
    legend.position = "inside",
    legend.position.inside = c(0.03, 0.96),
    legend.justification.inside = c(0, 1),
    legend.direction = "vertical",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 5.8),
    legend.spacing.y = unit(0.5, "pt"),
    plot.margin = margin(4, 5, 4, 5)
  ) +
  coord_cartesian(clip = "off")

output_dir <- file.path("examples", "minimal", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pdf_path <- file.path(output_dir, "synthetic_locus_demo.pdf")
png_path <- file.path(output_dir, "synthetic_locus_demo.png")

ggsave(
  filename = pdf_path,
  plot = plot_obj,
  width = 3.25,
  height = 1.95,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = png_path,
  plot = plot_obj,
  width = 3.25,
  height = 1.95,
  dpi = 300,
  bg = "white"
)

message("[DEMO] Saved compact synthetic panel to: ", png_path)
