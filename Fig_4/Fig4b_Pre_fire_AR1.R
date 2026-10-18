script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (is.null(script_file)) getwd() else dirname(normalizePath(script_file))
source(file.path(script_dir, "_plot_common.R"))

plot_shap_pair(
  csv_file = file.path(script_dir, "Fig_4_b_j.csv"),
  pair = 1,
  png_file = file.path(script_dir, "Fig4b_Pre_fire_AR1.png"),
  xlab = "Pre-fire AR(1)"
)
