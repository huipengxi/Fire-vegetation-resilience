# Figure 1b：纬度剖面

library(terra)
library(matrixStats)

trend <- rast("Sen_slopePerYear_annual.tif")

trend_mat <- as.matrix(trend, wide = TRUE)
n <- rowCounts(!is.na(trend_mat))

profile <- data.frame(
  Latitude = yFromRow(trend, seq_len(nrow(trend))),
  Mean = rowMeans2(trend_mat, na.rm = TRUE),
  SE = rowSds(trend_mat, na.rm = TRUE) / sqrt(n)
)

profile <- subset(profile, n > 0)
profile$CI_lower <- profile$Mean - 1.96 * profile$SE
profile$CI_upper <- profile$Mean + 1.96 * profile$SE

write.csv(
  profile[c("Latitude", "Mean", "CI_lower", "CI_upper")],
  "Figure1b_latitudinal_profile.csv",
  row.names = FALSE
)