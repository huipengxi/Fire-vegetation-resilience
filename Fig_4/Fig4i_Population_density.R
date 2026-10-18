script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (is.null(script_file)) getwd() else dirname(normalizePath(script_file))
source(file.path(script_dir, "_plot_common.R"))

plot_shap_pair(
  csv_file = file.path(script_dir, "Fig_4_b_j.csv"),
  pair = 8,
  png_file = file.path(script_dir, "Fig4i_Population_density.png"),
  xlab = expression(Population~density~(people~km^-2)),
  keep = function(x) x <= 400,
  x_limits = c(0, 400),
  x_breaks = seq(0, 400, by = 100),
  x_expand = c(0.02, 0.08)
)
