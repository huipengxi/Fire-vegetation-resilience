# Figure 1d: Latitudinal profile of the resilience baseline

library(terra)
library(matrixStats)

baseline <- rast("Resilience_baseline.tif")

baseline_mat <- as.matrix(baseline, wide = TRUE)
n <- rowCounts(!is.na(baseline_mat))

profile <- data.frame(
  Latitude = yFromRow(baseline, seq_len(nrow(baseline))),
  Mean = rowMeans2(baseline_mat, na.rm = TRUE),
  SE = rowSds(baseline_mat, na.rm = TRUE) / sqrt(n)
)

profile <- subset(profile, n > 1)

profile$CI_lower <- profile$Mean - 1.96 * profile$SE
profile$CI_upper <- profile$Mean + 1.96 * profile$SE

write.csv(
  profile[c("Latitude", "Mean", "CI_lower", "CI_upper")],
  "Figure1d_latitudinal_profile.csv",
  row.names = FALSE
)