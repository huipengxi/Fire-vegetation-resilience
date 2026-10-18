plot_shap_pair <- function(
    csv_file,
    pair,
    png_file,
    xlab,
    keep = function(x) rep(TRUE, length(x)),
    x_limits = NULL,
    x_breaks = ggplot2::waiver(),
    x_expand = c(0.02, 0.02),
    nbins = 50,
    shap_threshold = 1e-3
) {
  if (!file.exists(csv_file)) stop("File not found: ", csv_file)
  if (!pair %in% 1:9) stop("pair must be an integer from 1 to 9.")

  source_data <- read.csv(
    csv_file,
    skip = 1,
    check.names = FALSE,
    na.strings = c("", "NA")
  )

  index <- 2 * pair - 1
  df <- data.frame(
    Feature = source_data[[index]],
    SHAP = source_data[[index + 1]]
  )
  df <- df[
    is.finite(df$Feature) & is.finite(df$SHAP) & keep(df$Feature),
  ]

  df$bin <- cut(df$Feature, breaks = nbins, include.lowest = TRUE)
  groups <- split(df, df$bin, drop = TRUE)
  bin_summary <- do.call(
    rbind,
    lapply(groups, function(x) {
      data.frame(
        n_bin = nrow(x),
        n_explained = sum(abs(x$SHAP) > shap_threshold),
        x_median = median(x$Feature),
        y_median = median(x$SHAP)
      )
    })
  )
  bin_summary <- bin_summary[
    bin_summary$n_bin > 50 & bin_summary$n_explained >= 100,
  ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Feature, y = SHAP)) +
    ggplot2::geom_point(alpha = 0.2, size = 1, color = "black") +
    ggplot2::geom_line(
      data = bin_summary,
      ggplot2::aes(x = x_median, y = y_median),
      inherit.aes = FALSE,
      color = "red",
      linewidth = 1.2
    ) +
    ggplot2::scale_x_continuous(
      name = xlab,
      limits = x_limits,
      breaks = x_breaks,
      expand = ggplot2::expansion(mult = x_expand)
    ) +
    ggplot2::scale_y_continuous(
      name = "SHAP value [-]",
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 24, base_family = "Times New Roman") +
    ggplot2::theme(
      aspect.ratio = 1,
      text = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(size = 26, color = "black"),
      axis.title = ggplot2::element_text(size = 26, color = "black"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      axis.ticks.length = grid::unit(-0.3, "cm"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 1),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "black", fill = NA, linewidth = 1
      ),
      plot.margin = ggplot2::margin(10, 20, 10, 10, unit = "pt"),
      legend.position = "none"
    )

  ggplot2::ggsave(
    filename = png_file,
    plot = p,
    width = 15,
    height = 15,
    units = "cm",
    dpi = 600,
    bg = "white"
  )

  message(
    basename(png_file), ": displayed ", nrow(df),
    " rows; red-line bins ", nrow(bin_summary)
  )
  invisible(p)
}
