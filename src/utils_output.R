suppressPackageStartupMessages({
  library(glue)
  library(ggplot2)
})

save_plot_with_fallback <- function(plot_obj, pdf_path, png_path,
                                    plot_width, plot_height,
                                    png_scale = 0.8, png_dpi = 300) {
  if (file.exists(pdf_path)) file.remove(pdf_path)
  if (file.exists(png_path)) file.remove(png_path)

  pdf_saved <- tryCatch({
    ggsave(
      filename = pdf_path,
      plot = plot_obj,
      width = plot_width,
      height = plot_height,
      device = cairo_pdf,
      bg = "white"
    )
    TRUE
  }, error = function(e) {
    message(glue("[PLOT] PDF save failed: {basename(pdf_path)} -> {e$message}"))
    FALSE
  })

  if (!pdf_saved && file.exists(pdf_path)) {
    file.remove(pdf_path)
  }

  if (!pdf_saved) {
    ggsave(
      filename = png_path,
      plot = plot_obj,
      width = plot_width * png_scale,
      height = plot_height * png_scale,
      dpi = png_dpi,
      bg = "white"
    )
    return(list(format = "png", path = png_path))
  }

  list(format = "pdf", path = pdf_path)
}
