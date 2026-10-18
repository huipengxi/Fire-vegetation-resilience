library(ggplot2)

# Input and output files
data_file <- "Fig3c_Source_Data.csv"
figure_file <- "Fig3c_resilience_boxplot.png"

# Category order
area_levels <- c(
  "0–20%",
  "20–40%",
  "40–60%",
  "60–80%",
  "80–100%",
  ">100%"
)

frequency_levels <- c(
  "≤1 Fire",
  "1–3 Fires",
  "3–5 Fires",
  ">5 Fires"
)

# Fill colors
fire_colors <- c(
  "≤1 Fire" = "#4dac26",
  "1–3 Fires" = "#b8e186",
  "3–5 Fires" = "#f1b6da",
  ">5 Fires" = "#d01c8b"
)

# Read source data
plot_data <- read.csv(
  data_file,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8",
  check.names = FALSE
)

# Check required columns
required_columns <- c(
  "Burned_area_fraction",
  "Fire_frequency",
  "Mean_change_pct",
  "Lower_whisker_pct",
  "Q1_change_pct",
  "Median_change_pct",
  "Q3_change_pct",
  "Upper_whisker_pct"
)

missing_columns <- setdiff(required_columns, names(plot_data))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Set category order
plot_data$Burned_area_fraction <- factor(
  plot_data$Burned_area_fraction,
  levels = area_levels
)

plot_data$Fire_frequency <- factor(
  plot_data$Fire_frequency,
  levels = frequency_levels
)

# Remove undefined categories
plot_data <- plot_data[
  !is.na(plot_data$Burned_area_fraction) &
    !is.na(plot_data$Fire_frequency),
]

# Align boxes and mean points
dodge_position <- position_dodge(width = 0.8)

# Draw the figure
p <- ggplot(
  plot_data,
  aes(
    x = Burned_area_fraction,
    fill = Fire_frequency,
    group = interaction(Burned_area_fraction, Fire_frequency)
  )
) +
  geom_boxplot(
    aes(
      ymin = Lower_whisker_pct,
      lower = Q1_change_pct,
      middle = Median_change_pct,
      upper = Q3_change_pct,
      ymax = Upper_whisker_pct
    ),
    stat = "identity",
    position = dodge_position,
    width = 0.65,
    color = "black",
    linewidth = 0.5,
    alpha = 1
  ) +
  geom_point(
    aes(y = Mean_change_pct),
    position = dodge_position,
    shape = 21,
    size = 2.5,
    stroke = 0.6,
    fill = "white",
    color = "black"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = fire_colors,
    drop = FALSE
  ) +
  scale_x_discrete(
    limits = area_levels,
    drop = TRUE,
    expand = expansion(add = c(0.55, 0.55))
  ) +
  scale_y_continuous(
    breaks = c(-100, 0, 100),
    minor_breaks = c(-150, -50, 50, 150),
    expand = expansion(mult = c(0.02, 0.04)),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  coord_cartesian(
    ylim = c(-150, 150),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = "Change in resilience (%)",
    fill = NULL
  ) +
  theme_minimal(
    base_size = 24,
    base_family = "Times New Roman"
  ) +
  theme(
    text = element_text(
      family = "Times New Roman",
      color = "black"
    ),
    axis.text.x = element_text(
      size = 26,
      angle = 0,
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 26,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 26,
      color = "black",
      margin = margin(r = 10)
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 1
    ),
    axis.ticks.length.x.bottom = grid::unit(-0.45, "cm"),
    axis.ticks.length.y.left = grid::unit(-0.40, "cm"),
    axis.minor.ticks.y.left = element_line(
      color = "black",
      linewidth = 0.7
    ),
    axis.minor.ticks.length.y.left = rel(0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "gray85",
      linetype = "dotted",
      linewidth = 0.5
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    panel.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
  )

# Display and save the figure
print(p)

ggsave(
  filename = figure_file,
  plot = p,
  width = 1500,
  height = 550,
  units = "px",
  dpi = 300,
  bg = "white"
)
