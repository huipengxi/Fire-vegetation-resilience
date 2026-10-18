suppressPackageStartupMessages({
  library(ggplot2)
  library(showtext)
})


dat <- read.csv(
  file.path("Fig2c_Source_Data.csv"),
  check.names = FALSE
)

dat$Vegetation_type <- factor(
  dat$Vegetation_type,
  levels = c("BF", "MF", "Shrub", "Sav", "Grass")
)
dat$Time <- factor(
  dat$Time,
  levels = c(
    "Pre-fire",
    "1 year post-fire",
    "3 years post-fire",
    "5 years post-fire"
  )
)
dat <- dat[order(dat$Vegetation_type, dat$Time), ]

# Read the original fill colours directly from Source Data.
time_colours <- setNames(dat$Fill_colour, dat$Time)
time_colours <- time_colours[!duplicated(names(time_colours))]

font_add(
  family = "Times New Roman",
  regular = "C:/Windows/Fonts/times.ttf",
  bold = "C:/Windows/Fonts/timesbd.ttf",
  italic = "C:/Windows/Fonts/timesi.ttf",
  bolditalic = "C:/Windows/Fonts/timesbi.ttf"
)
showtext_auto()

dodge <- position_dodge(width = 0.80)

p <- ggplot(
  dat,
  aes(
    x = Vegetation_type,
    fill = Time,
    group = interaction(Vegetation_type, Time)
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
    position = dodge,
    width = 0.65,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_point(
    aes(y = Mean_change_pct),
    position = dodge,
    shape = 21,
    size = 2.5,
    stroke = 0.6,
    fill = "white",
    colour = "black"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = time_colours) +
  scale_y_continuous(
    breaks = c(-100, 0, 100),
    minor_breaks = c(-150, -50, 50, 150),
    expand = expansion(mult = c(0.02, 0.04)),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  coord_cartesian(ylim = c(-185, 185), clip = "off") +
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
      colour = "black"
    ),
    axis.text.x = element_text(size = 26, colour = "black"),
    axis.text.y = element_text(size = 26, colour = "black"),
    axis.title.y = element_text(
      size = 26,
      colour = "black",
      margin = margin(r = 10)
    ),
    axis.ticks = element_line(colour = "black", linewidth = 1),
    axis.ticks.length.x.bottom = grid::unit(-0.45, "cm"),
    axis.ticks.length.y.left = grid::unit(-0.40, "cm"),
    axis.minor.ticks.y.left = element_line(
      colour = "black",
      linewidth = 0.7
    ),
    axis.minor.ticks.length.y.left = rel(0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      colour = "grey85",
      linetype = "dotted",
      linewidth = 0.5
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    ),
    legend.position = "none",
    plot.margin = margin(10, 15, 10, 10)
  )

print(p)

ggsave(
  file.path(data_dir, "Fig2c_boxplot_1500x500.png"),
  p,
  width = 1500,
  height = 500,
  units = "px",
  dpi = 300,
  device = ragg::agg_png,
  bg = "white"
)

ggsave(
  file.path(data_dir, "Fig2c_boxplot.pdf"),
  p,
  width = 15,
  height = 5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
