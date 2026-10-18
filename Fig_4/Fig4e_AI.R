script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (is.null(script_file)) getwd() else dirname(normalizePath(script_file))
source(file.path(script_dir, "_plot_common.R"))

plot_shap_pair(
  csv_file = file.path(script_dir, "Fig_4_b_j.csv"),
  pair = 4,
  png_file = file.path(script_dir, "Fig4e_AI.png"),
  xlab = "Aridity index",
  keep = function(x) x <= 2
)
