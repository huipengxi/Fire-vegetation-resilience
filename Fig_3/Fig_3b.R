library(ggplot2)
library(showtext)
library(grid)

# Read the final source data
plot_data <- read.csv(
  paste0(
    "D:/R-code/P2/20261006/SourceData/Figure_3/b/",
    "Source_Data_Fig3b.csv"
  ),
  check.names = FALSE
)

# Recover the common coordinate ranges
fraction_min <- plot_data$Fraction_axis_min[1]
fraction_max <- plot_data$Fraction_axis_max[1]
frequency_min <- plot_data$Frequency_axis_min[1]
frequency_max <- plot_data$Frequency_axis_max[1]

# Fixed integer axis labels
fraction_breaks_original <- c(0, 60, 120)
frequency_breaks_original <- c(1, 3, 5, 7, 9)

fraction_breaks_norm <- (
  fraction_breaks_original - fraction_min
) / (
  fraction_max - fraction_min
)

frequency_breaks_norm <- (
  frequency_breaks_original - frequency_min
) / (
  frequency_max - frequency_min
)

# Set the font
font_add(
  family = "Times New Roman",
  regular = "C:/Windows/Fonts/times.ttf"
)
showtext_auto(enable = TRUE)

# Shared color palette
my_colors <- c(
  "#b35806",
  "#f1a340",
  "#fee0b6",
  "#d8daeb",
  "#998ec3",
  "#542788"
)

panel_padding <- 0.04

# Draw Fig. 3b
p <- ggplot(plot_data) +
  geom_rect(
    aes(
      xmin = Fraction_norm_lower,
      xmax = Fraction_norm_upper,
      ymin = Frequency_norm_lower,
      ymax = Frequency_norm_upper,
      fill = Displayed_Maximum_Recovery_Time_months
    ),
    color = NA
  ) +
  scale_fill_gradientn(
    colors = my_colors,
    name = "Maximum recovery time\n(months)",
    limits = c(0, 36),
    oob = scales::squish,
    breaks = c(0, 10, 20, 30),
    labels = c("0", "10", "20", "30"),
    guide = guide_colorbar(
      title.position = "right",
      title.hjust = 0.5
    )
  ) +
  scale_x_continuous(
    name = "Mean annual cumulative\nburned fraction (%)",
    breaks = fraction_breaks_norm,
    labels = sprintf("%.0f", fraction_breaks_original),
    limits = c(-panel_padding, 1 + panel_padding),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Mean annual burn frequency",
    breaks = frequency_breaks_norm,
    labels = sprintf("%.0f", frequency_breaks_original),
    limits = c(-panel_padding, 1 + panel_padding),
    expand = c(0, 0)
  ) +
  coord_fixed(ratio = 1, clip = "on") +
  theme_minimal(base_size = 24) +
  theme(
    text = element_text(
      family = "Times New Roman",
      color = "black"
    ),
    axis.title.x = element_text(
      size = 26,
      margin = margin(t = 14)
    ),
    axis.title.y = element_text(
      size = 26,
      margin = margin(r = 14)
    ),
    axis.text.x = element_text(
      size = 26,
      color = "black",
      margin = margin(t = 10)
    ),
    axis.text.y = element_text(
      size = 26,
      color = "black",
      margin = margin(r = 10)
    ),
    axis.ticks.length = unit(-0.30, "cm"),
    axis.ticks = element_line(
      color = "black",
      linewidth = 1
    ),
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    legend.title = element_text(
      size = 24,
      angle = 270,
      vjust = 0.5
    ),
    legend.text = element_text(
      size = 20,
      color = "black"
    ),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.6, "cm")
  )

p
