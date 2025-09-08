#####################################################################################################################
##
## This code was used to conduct analysis for the paper "Fire-driven declines in vegetation resilience across the tropics".
##
#####################################################################################################################

library(stats)
library(sp)
library(raster)
library(terra)
library(trend)
library(xgboost)
library(SHAPforxgboost)
library(SHAPforxgboost)
library(matrixStats)
library(ggplot2)
library(dplyr)
library(car)

#********************************************
#-----1.Read the raw data and detrend it-----
#********************************************

Data <- read.csv("./ExampleData.csv")
Data_TS <- Data[,2]
ts_Data_STL <- ts(Data_TS, start = c(2001, 1), frequency = 12)
decomposed <- stl(ts_Data_STL, s.window = "periodic", t.window = 23, l.window = 25, robust = FALSE)

seasonal <- decomposed$time.series[, "seasonal"] # extract the seasonal trend
seasonal_kNDVI <- data.frame(seasonal)
trend <- decomposed$time.series[, "trend"] # extract the long-term trend
trend_kNDVI <- data.frame(trend)

remainder <- decomposed$time.series[, "remainder"] # extract the residuals
plot(decomposed)

#**********************************************************************************
#-----2.Calculate the autocorrelation coefficient AR(1) based on the residuals-----
#**********************************************************************************

window_size <- 60
calc_ar1_for_pixel <- function(pixel_ts, window_size) {
  if (length(pixel_ts) < window_size) {
    return(NA)
  }
  ar1_values <- sapply(1:(length(pixel_ts) - window_size + 1), function(t) {
    window_data <- pixel_ts[t:(t + window_size - 1)]
    ar_model <- ar.ols(window_data, aic = FALSE, order.max = 1)
    return(ar_model$ar)
  })
  return(ar1_values)
}

pixel_ts <- remainder
# Calculate the AR(1)
ar1_values <- calc_ar1_for_pixel(pixel_ts, window_size)

#****************************************************************************************************************
#-----3. Trend analysis was conducted using the Mann¨CKendall test in combination with Sen¡¯s slope estimator.-----
#****************************************************************************************************************

fl <- list.files(path = "./AR1_kNDVI/5_yr/Yearly",
                 pattern = ".tif$",full.names=TRUE)

firs <- rast(fl)
fun_sen <- function(x){
  if(length(na.omit(x)) < 17)
    return(c(NA,NA,NA))
  MK_estimate <- trend::sens.slope(ts(na.omit(x), start=2006, end=2022, frequency=1), conf.level=0.95)
  slope <- MK_estimate$estimate
  MK_test <- MK_estimate$p.value
  Zs <- MK_estimate$statistic
  return(c(Zs,slope,MK_test))
}

firs_sen <- app(firs,fun_sen,cores=8)
names(firs_sen) <- c("Z","slope","p-value")

plot(firs_sen)
writeRaster(firs_sen,filename="./firs_sen_kNDVI.tif",names=firs_sen@ptr[["names"]])

# Extracting the latitude profiles of resilience trends

r <- raster("./Slope_Reslience_GIMMS_kNDVI.tif")
x_mat <- as.matrix(r)

mean_vals <- rowMeans2(x_mat, na.rm = TRUE)
sd_vals   <- rowSds(x_mat, na.rm = TRUE)
n_vals    <- rowCounts(!is.na(x_mat))
se_vals   <- sd_vals / sqrt(n_vals)

ci_upper <- mean_vals + 1.96 * se_vals
ci_lower <- mean_vals - 1.96 * se_vals

lat_start <- 30
lat_step  <- -0.05
lat_vals  <- seq(lat_start, lat_start + (nrow(x_mat)-1) * lat_step, by = lat_step)

df_lat <- data.frame(
  Latitude = lat_vals,
  Mean     = mean_vals,
  CI_lower = ci_lower,
  CI_upper = ci_upper
)
write.csv(df_lat, "./Lat_plot/GIMMS_kNDVI/Lon_trend_stat.csv", row.names = FALSE)


#*******************************************************
#-----4. Integrate all dynamic and static variables-----
#*******************************************************

#*****¦¤Resilience*****
pre_fire_year <- 3
years <- 2009:2020
Res_list <- list()
for (yr in years) {
  out_dir <- file.path("./AR1_Diff/", paste0("post_", pre_fire_year))
  file_path <- file.path(out_dir, paste0("AR1_Diff_", yr, "_masked.tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  Res_Diff <- rast(file_path)
  Res_list <- c(Res_list, Res_Diff)
}

Res_Diff_stack <- rast(Res_list)

names(Res_Diff_stack) <- paste0("AR1_Diff_", 2009:2020)

#*****¦¤Albedo*****
years <- 2009:2020
Alb_list <- list()
for (yr in years) {
  out_dir <- file.path("./Albedo_change/")
  file_path <- file.path(out_dir, paste0("merged_diff_postfire_Alb_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  Alb_Diff <- rast(file_path)
  Alb_list <- c(Alb_list, Alb_Diff)
}
Alb_Diff_stack <- rast(Alb_list)

#*****¦¤ET*****
years <- 2009:2020
ET_list <- list()
for (yr in years) {
  out_dir <- file.path("./ET_change/")
  file_path <- file.path(out_dir, paste0("merged_diff_postfire_ET_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  ET_Diff <- rast(file_path)
  ET_list <- c(ET_list, ET_Diff)
}

ET_Diff_stack <- rast(ET_list)

#*****¦¤PRE*****
years <- 2009:2020
PRE_list <- list()
for (yr in years) {
  out_dir <- file.path("./PRE_change/")
  file_path <- file.path(out_dir, paste0("merged_diff_postfire_PRE_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  PRE_Diff <- rast(file_path)
  PRE_list <- c(PRE_list, PRE_Diff)
}
PRE_Diff_stack <- rast(PRE_list)


#*****¦¤SM*****
years <- 2009:2020
SM_list <- list()
for (yr in years) {
  out_dir <- file.path("./SM_change/")
  file_path <- file.path(out_dir, paste0("merged_diff_postfire_SM_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  SM_Diff <- rast(file_path)
  SM_list <- c(SM_list, SM_Diff)
}
SM_Diff_stack <- rast(SM_list)


#*****¦¤SW*****
years <- 2009:2020
SW_list <- list()
for (yr in years) {
  out_dir <- file.path("./SW_change/merged/")
  file_path <- file.path(out_dir, paste0("merged_diff_postfire_SW_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  SW_Diff <- rast(file_path)
  SW_list <- c(SW_list, SW_Diff)
}
SW_Diff_stack <- rast(SW_list)

#*****¦¤TMP*****
years <- 2009:2020
TMP_list <- list()
for (yr in years) {
  out_dir <- file.path("./TMP_change/")
  file_path <- file.path(out_dir, paste0("merged_diff_postfire_TMP_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  TMP_Diff <- rast(file_path)
  TMP_list <- c(TMP_list, TMP_Diff)
}
TMP_Diff_stack <- rast(TMP_list)

#*****¦¤VPD*****
years <- 2009:2020
VPD_list <- list()
for (yr in years) {
  out_dir <- file.path("./VPD_change/")
  file_path <- file.path(out_dir, paste0("merged_diff_postfire_VPD_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  VPD_Diff <- rast(file_path)
  VPD_list <- c(VPD_list, VPD_Diff)
}
VPD_Diff_stack <- rast(VPD_list)

#*****dNBR*****
years <- 2009:2020
dNBR_list <- list()
for (yr in years) {
  out_dir <- file.path("./Burn_dNBR/Mask/")
  file_path <- file.path(out_dir, paste0("dNBR_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  dNBR_Diff <- rast(file_path)
  dNBR_list <- c(dNBR_list, dNBR_Diff)
}
dNBR_Diff_stack <- rast(dNBR_list)


#*****dLAI*****
years <- 2009:2020
dLAI_list <- list()
for (yr in years) {
  out_dir <- file.path("./Burn_dLAI/Mask/")
  file_path <- file.path(out_dir, paste0("dLAI_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  dLAI_Diff <- rast(file_path)
  dLAI_list <- c(dLAI_list, dLAI_Diff)
}
dLAI_Diff_stack <- rast(dLAI_list)


#*****Pre-fire AR(1)*****
years <- 2009:2020
Pre_Res_list <- list()
for (yr in years) {
  out_dir <- file.path("./Pre_3_AR(1)/")
  file_path <- file.path(out_dir, paste0("AR1_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  Pre_Res_Diff <- rast(file_path)
  Pre_Res_list <- c(Pre_Res_list, Pre_Res_Diff)
}
Pre_Res_Diff_stack <- rast(Pre_Res_list)

#*****Pre-fire ¦¤kNDVI*****
# The percent of rainy seasons greater than baseline in the 36 months before the fire

years <- 2009:2020
Per_Res_list <- list()
for (yr in years) {
  out_dir <- file.path("./percent_kNDVI/")
  file_path <- file.path(out_dir, paste0("Percent_kNDVI_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  Pre_Res_Diff <- rast(file_path)
  Per_Res_list <- c(Per_Res_list, Pre_Res_Diff)
}
Percent_Diff_stack <- rast(Per_Res_list)

#*****Footprint*****
years <- 2009:2020
Pop_list <- list()
for (yr in years) {
  out_dir <- file.path("./Footprint/")
  file_path <- file.path(out_dir, paste0("Pop_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  Pop_Diff <- rast(file_path)
  Pop_list <- c(Pop_list, Pop_Diff)
}
Pop_stack <- rast(Pop_list)


#*****SOC*****
SOC <- rast("./Soils/SOC/SOC_content.tif")

#*****Soil Clay*****
Clay <- rast("./Soils/Clay/Clay_content.tif")

#*****Depth to bedrock*****
Soil_Depth <- rast("./Soils/soil_depth.tif")

#*****AI*****
AI <- rast("./AI_Data/AI_averge_5km.tif")

#*****AGB*****
AGB_2010 <- rast("./Benchmark_AGB_2010.tif")

#*****Burn area*****
years <- 2009:2020
Burn_area_list <- list()
for (yr in years) {
  out_dir <- file.path("./Burn_Area/")
  file_path <- file.path(out_dir, paste0("Burn_Area_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  Burn_area_Diff <- rast(file_path)
  Burn_area_list <- c(Burn_area_list, Burn_area_Diff)
}
Burn_area_Diff_stack <- rast(Burn_area_list)

##*****Burn frequency*****
years <- 2009:2020
Burn_freq_list <- list()
for (yr in years) {
  out_dir <- file.path("./Fire_Frequency/")
  file_path <- file.path(out_dir, paste0("Burn_Frequency_", yr, ".tif"))
  if (!file.exists(file_path)) {
    warning("Missing file£º", file_path)
    next
  }
  Burn_freq_Diff <- rast(file_path)
  Burn_freq_list <- c(Burn_freq_list, Burn_freq_Diff)
}
Burn_freq_Diff_stack <- rast(Burn_freq_list)



# === Dynamic variables ===
diff_stacks <- list(
  Alb      = Alb_Diff_stack,
  Res      = Res_Diff_stack,
  ET       = ET_Diff_stack,
  PRE      = PRE_Diff_stack,
  SM       = SM_Diff_stack,
  SW       = SW_Diff_stack,
  TMP      = TMP_Diff_stack,
  VPD      = VPD_Diff_stack,
  dNBR     = dNBR_Diff_stack,
  dLAI     = dLAI_Diff_stack,
  Pre_Res = Pre_Res_Diff_stack,
  Per_Res  = Pre_Res_Diff_stack,
  Burn_area = Burn_area_Diff_stack,
  Burn_freq = Burn_freq_Diff_stack,
  Pop_density = Pop_stack
)

# === Static variables ===
static_rasters <- list(
  SOC         = SOC,
  AI          = AI,
  Clay        = Clay,
  Soil_Depth  = Soil_Depth,
  AGB         = AGB
)

years <- 2009:2020
n_years <- length(years)
df_list <- list()

for (i in seq_along(years)) {
  yr <- years[i]
  message("Processing year: ", yr)
  
  n_pixel <- ncell(SOC)
  vars_df <- data.frame(Year = rep(yr, n_pixel))
  
  # === Dynamic variables ===
  for (var_name in names(diff_stacks)) {
    layer_i <- diff_stacks[[var_name]][[i]]
    val_i <- values(layer_i)
    
    if (length(val_i) != n_pixel) {
      warning("jump", var_name, " in ", yr, " year£ºpixels are inconsistent")
      val_i <- rep(NA, n_pixel)
    }
    vars_df[[var_name]] <- val_i
  }
  
  # === Static variables ===
  for (sname in names(static_rasters)) {
    val_static <- values(static_rasters[[sname]])
    if (length(val_static) != n_pixel) {
      warning("Static variables ", sname, "size abnormal")
      val_static <- rep(NA, n_pixel)
    }
    vars_df[[sname]] <- val_static 
  }
  vars_df <- vars_df %>% filter(if_all(everything(), ~ !is.na(.)))
  df_list[[i]] <- vars_df
}

final_df <- bind_rows(df_list)
print(head(final_df))
write.csv(final_df, "./final_diff_variables_2009_2020.csv", row.names = FALSE)

#**************************************************************************
#-----5. XGBoost model for identifying critical environmental controls-----
#**************************************************************************

rm(list = ls())
gc()

final_df <- read.csv("./final_diff_variables_2009_2020.csv")
df_no_year <- final_df[, !names(final_df) %in% "Year"]
vif_model <- lm(Res ~ ., data = df_no_year)
vif_vals <- vif(vif_model)
print(vif_vals)
high_vif <- names(vif_vals[vif_vals > 3])
cat("Variables with high collinearity (VIF > 3):\n")
print(high_vif) # The results show that the VIF of all variables is less than 3

repeat_xgb_cv <- function(data, 
                          n_repeats = 50, 
                          nfold = 3, 
                          sample_size = 10000,
                          nrounds = 1000,
                          early_stop = 30,
                          seed = 123,
                          replace = TRUE) {
  set.seed(seed)
  
  r2_train_vec <- numeric(n_repeats)
  r2_test_vec <- numeric(n_repeats)
  shap_list <- list()
  
  for (i in 1:n_repeats) {
    cat("Running iteration", i, "...\n")
    train_idx <- sample(seq_len(nrow(data)), size = sample_size)
    train_data <- data[train_idx, ]
    test_data <- data[-train_idx, ]
    
    train_matrix <- xgb.DMatrix(data = as.matrix(train_data[, !names(train_data) %in% "Res"]),
                                label = train_data$Res)
    test_matrix <- xgb.DMatrix(data = as.matrix(test_data[, !names(test_data) %in% "Res"]),
                               label = test_data$Res)
    
    cv_model <- xgb.cv(
      params = list(
        objective = "reg:squarederror",
        eval_metric = "rmse",
        eta = 0.03,
        max_depth = 10,
        subsample = 0.7,
        colsample_bytree = 0.7,
        min_child_weight = 10,
        gamma = 0.3,
        alpha = 5,
        lambda = 3
      ),
      data = train_matrix,
      nrounds = nrounds,
      nfold = nfold,
      early_stopping_rounds = early_stop,
      verbose = 0
    )
    
    best_n <- cv_model$best_iteration
    
    model <- xgb.train(
      params = cv_model$params,
      data = train_matrix,
      nrounds = best_n,
      verbose = 0
    )
    
    pred_train <- predict(model, train_matrix)
    pred_test  <- predict(model, test_matrix)
    
    r2_train <- 1 - sum((train_data$Res - pred_train)^2) / sum((train_data$Res - mean(train_data$Res))^2)
    r2_test  <- 1 - sum((test_data$Res - pred_test)^2) / sum((test_data$Res - mean(test_data$Res))^2)
    
    r2_train_vec[i] <- r2_train
    r2_test_vec[i]  <- r2_test
    
    # SHAP
    shap <- shap.values(xgb_model = model,
                        X_train = as.matrix(train_data[, !names(train_data) %in% "Res"]))$mean_shap_score
    shap_list[[i]] <- shap
  }
  
  shap_mat <- do.call(cbind, shap_list)
  shap_mean <- rowMeans(shap_mat)
  shap_df <- data.frame(Variable = names(shap_mean), Mean_SHAP = shap_mean)
  shap_df <- shap_df[order(-shap_df$Mean_SHAP), ]
  
  return(list(
    r2_train = r2_train_vec,
    r2_test = r2_test_vec,
    shap_df = shap_df,
    shap_matrix = shap_mat
  ))
}
results <- repeat_xgb_cv(data = final_df, n_repeats = 100)
saveRDS(results, file = "./results_xgb_cv100_2.rds")
