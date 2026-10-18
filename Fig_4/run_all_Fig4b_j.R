script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (is.null(script_file)) getwd() else dirname(normalizePath(script_file))

scripts <- c(
  "Fig4b_Pre_fire_AR1.R",
  "Fig4c_AGB.R",
  "Fig4d_ET.R",
  "Fig4e_AI.R",
  "Fig4f_Clay.R",
  "Fig4g_VPD.R",
  "Fig4h_Soil_depth.R",
  "Fig4i_Population_density.R",
  "Fig4j_dLAI.R"
)

for (script in scripts) {
  message("Running: ", script)
  source(file.path(script_dir, script), local = new.env(parent = globalenv()))
}

message("All nine figures saved in: ", script_dir)
