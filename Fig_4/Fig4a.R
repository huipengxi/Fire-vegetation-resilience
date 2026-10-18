suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(showtext)
})

input_file <- "Fig4a.csv"
output_file <- file.path(dirname(input_file), "Fig4a.png")

font_add(
  family = "Times New Roman",
  regular = "C:/Windows/Fonts/times.ttf"
)
showtext_auto()
showtext_opts(dpi = 600)

# Source Data: 19 variables x 100 model runs.
# The current CSV is saved in Windows Chinese encoding (GBK/GB18030).
csv_text <- iconv(
  readLines(input_file, encoding = "bytes"),
  from = "GB18030",
  to = "UTF-8"
)
shap_data <- read.csv(
  text = paste(csv_text, collapse = "\n"),
  check.names = FALSE
)

required_columns <- c("Run", "Variable", "SHAP_pct")
if (!all(required_columns %in% names(shap_data))) {
  stop("Fig4a.csv must contain: Run, Variable, and SHAP_pct.")
}

shap_data <- shap_data %>%
  select(all_of(required_columns)) %>%
  filter(is.finite(SHAP_pct))

# Bar height = mean; error bar = +/- SD across the 100 runs.
shap_summary <- shap_data %>%
  group_by(Variable) %>%
  summarise(
    Mean_SHAP = mean(SHAP_pct),
    SD_SHAP = sd(SHAP_pct),
    .groups = "drop"
  ) %>%
  arrange(Mean_SHAP)

# Use the mean contribution to determine the bottom-to-top order.
variable_levels <- shap_summary$Variable
shap_summary$Variable <- factor(shap_summary$Variable, levels = variable_levels)
shap_data$Variable <- factor(shap_data$Variable, levels = variable_levels)

colour_values <- setNames(
  colorRampPalette(c("#a1d99b", "#00441b"))(nrow(shap_summary)),
  variable_levels
)

r2_value <- 0.64

p <- ggplot(
  shap_summary,
  aes(x = Variable, y = Mean_SHAP, fill = Variable)
) +
  geom_col(width = 0.88) +
  geom_point(
    data = shap_data,
    aes(x = Variable, y = SHAP_pct),
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0.15,
      height = 0,
      seed = 20261009
    ),
    shape = 21,
    size = 1.15,
    stroke = 0.18,
    colour = "#202020",
    fill = "white",
    alpha = 0.55
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_SHAP - SD_SHAP,
      ymax = Mean_SHAP + SD_SHAP
    ),
    width = 0.30,
    linewidth = 0.50
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = colour_values) +
  scale_y_continuous(
    name = "SHAP percentage contribution (%)",
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  labs(x = NULL) +
  annotate(
    "text",
    x = 1.5,
    y = max(shap_summary$Mean_SHAP) * 0.55,
    label = sprintf("R\u00B2 = %.2f", r2_value),
    family = "Times New Roman",
    size = 8,
    hjust = 0
  ) +
  theme_minimal(
    base_size = 22,
    base_family = "Times New Roman"
  ) +
  theme(
    text = element_text(
      family = "Times New Roman",
      colour = "black"
    ),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.text.x = element_text(size = 22, colour = "black"),
    axis.title.x = element_text(
      size = 24,
      margin = margin(t = 10),
      colour = "black"
    ),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    ),
    axis.ticks = element_line(colour = "black", linewidth = 1),
    axis.ticks.length = grid::unit(-0.2, "cm"),
    legend.position = "none"
  )

print(p)

ggsave(
  filename = output_file,
  plot = p,
  width = 7,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

message("Saved: ", output_file)
