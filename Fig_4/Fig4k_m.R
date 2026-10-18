suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

input_csv <- "Fig4k-m_Source_Data.csv"
output_dir <- dirname(input_csv)

dat <- read.csv(input_csv, check.names = FALSE, fileEncoding = "UTF-8")

required_columns <- c(
  "Figure_panel", "Predictor", "Quantile_tau",
  "Standardized_coefficient", "Bootstrap_SE", "t_statistic",
  "Residual_df", "P_value", "CI_95_lower", "CI_95_upper", "Significance"
)
if (!all(required_columns %in% names(dat))) {
  stop("The Source Data file does not contain all required columns.")
}

predictor_order <- c(
  "AGB", "Pre-fire AR(1)",
  "Delta ET", "Delta ET squared", "Delta VPD", "Delta VPD squared",
  "Soil clay", "Depth to bedrock"
)
dat$Predictor <- factor(dat$Predictor, levels = predictor_order)

colours <- c(
  "AGB" = "#D01C8B",
  "Pre-fire AR(1)" = "#0571B0",
  "Delta ET" = "#E66101",
  "Delta ET squared" = "#FDB863",
  "Delta VPD" = "#5E3C99",
  "Delta VPD squared" = "#B8AEDB",
  "Soil clay" = "#4CAF50",
  "Depth to bedrock" = "#7F8C8D"
)

legend_labels <- c(
  "AGB" = "AGB",
  "Pre-fire AR(1)" = "Pre-fire AR(1)",
  "Delta ET" = "ET",
  "Delta ET squared" = "ET^2",
  "Delta VPD" = "VPD",
  "Delta VPD squared" = "VPD^2",
  "Soil clay" = "Soil clay",
  "Depth to bedrock" = "Depth to bedrock"
)

theme_qr <- function() {
  theme_classic(base_size = 11, base_family = "sans") +
    theme(
      axis.title = element_text(size = 11, colour = "black"),
      axis.text = element_text(size = 9, colour = "black"),
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.position = "bottom",
      legend.spacing.x = unit(2, "pt"),
      legend.key.width = unit(18, "pt"),
      plot.title = element_text(size = 12, face = "bold"),
      plot.tag = element_text(size = 13, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

make_panel <- function(panel, title, tag, legend_rows = 1) {
  panel_data <- filter(dat, Figure_panel == panel)

  ggplot(
    panel_data,
    aes(
      x = Quantile_tau,
      y = Standardized_coefficient,
      colour = Predictor,
      fill = Predictor,
      group = Predictor
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.4
    ) +
    geom_ribbon(
      aes(ymin = CI_95_lower, ymax = CI_95_upper),
      alpha = 0.12,
      colour = NA
    ) +
    geom_line(linewidth = 1) +
    geom_point(
      aes(shape = Significance),
      size = 2,
      stroke = 0.8,
      fill = "white"
    ) +
    scale_colour_manual(values = colours, labels = legend_labels, drop = TRUE) +
    scale_fill_manual(values = colours, labels = legend_labels, drop = TRUE) +
    scale_shape_manual(
      values = c("Significant" = 16, "Not significant" = 1)
    ) +
    scale_x_continuous(
      breaks = c(0.1, 0.3, 0.5, 0.7, 0.9),
      limits = c(0.1, 0.9)
    ) +
    labs(
      x = expression("Resilience quantile (" * tau * ")"),
      y = "Standardized coefficient",
      title = title,
      tag = tag
    ) +
    guides(
      shape = "none",
      colour = guide_legend(nrow = legend_rows, byrow = TRUE),
      fill = guide_legend(nrow = legend_rows, byrow = TRUE)
    ) +
    theme_qr()
}

p_k <- make_panel(
  "Fig. 4k",
  "Pre-fire conditions and biomass",
  "k"
)
p_l <- make_panel(
  "Fig. 4l",
  "Water-stress effects",
  "l",
  legend_rows = 2
)
p_m <- make_panel(
  "Fig. 4m",
  "Soil properties",
  "m",
  legend_rows = 2
)

draw_three_panels <- function() {
  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(layout = grid::grid.layout(1, 3))
  )
  print(p_k, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_l, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p_m, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
}

# Preview in the active graphics device.
draw_three_panels()

# PNG: 19 cm wide, 600 dpi.
png(
  file.path(output_dir, "Fig4k-m_quantile_regression.png"),
  width = 19,
  height = 8,
  units = "cm",
  res = 600,
  bg = "white"
)
draw_three_panels()
dev.off()

# Vector PDF using the same physical dimensions.
cairo_pdf(
  file.path(output_dir, "Fig4k-m_quantile_regression.pdf"),
  width = 19 / 2.54,
  height = 8 / 2.54,
  bg = "white"
)
draw_three_panels()
dev.off()

