# Main analysis code for “Declining resilience to fire disturbances across tropical ecosystems”

library(terra)
library(trend)
library(matrixStats)
library(xgboost)
library(SHAPforxgboost)
library(car)


# 1. Detrend the monthly kNDVI series
data <- read.csv("./ExampleData.csv")
kndvi_ts <- ts(data[[2]], start = c(2001, 1), frequency = 12)
kndvi_stl <- stl(
  kndvi_ts,
  s.window = "periodic",
  t.window = 23,
  l.window = 25
)

seasonal_kNDVI <- data.frame(seasonal = kndvi_stl$time.series[, "seasonal"])
trend_kNDVI <- data.frame(trend = kndvi_stl$time.series[, "trend"])
remainder <- kndvi_stl$time.series[, "remainder"]


# 2. Calculate AR(1) in a moving window
moving_ar1 <- function(x, window = 60) {
  if (length(x) < window) return(NA_real_)

  vapply(seq_len(length(x) - window + 1), function(i) {
    fit <- ar.ols(x[i:(i + window - 1)], aic = FALSE, order.max = 1)
    unname(fit$ar[1])
  }, numeric(1))
}

ar1_values <- moving_ar1(remainder)


# 3. XGBoost and SHAP analysis
df <- Model_data
df <- df[
  is.finite(df$Res) &
    df$Res > 0 &
    complete.cases(df),
]

response <- "Res"
predictors <- setdiff(names(df), response)

stopifnot(
  response %in% names(df),
  nrow(df) > 1,
  all(vapply(df[, predictors, drop = FALSE],
             is.numeric,
             logical(1)))
)

xgb_params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  eta = 0.03,
  max_depth = 10,
  subsample = 0.7,
  colsample_bytree = 0.7,
  min_child_weight = 10,
  gamma = 0.3,
  alpha = 5,
  lambda = 3,
  nthread = 1
)

repeat_xgb <- function(
    data,
    response,
    predictors,
    params,
    n_repeats = 100,
    nfold = 3,
    sample_size = 50000,
    nrounds = 1000,
    early_stopping_rounds = 30,
    seed = 123
) {
  
  n <- nrow(data)
  metrics <- data.frame(
    Repeat = seq_len(n_repeats),
    Best_iteration = NA_integer_,
    Train_R2 = NA_real_,
    Test_R2 = NA_real_
  )
  
  shap_list <- vector("list", n_repeats)
  
  calculate_r2 <- function(observed, predicted) {
    1 - sum((observed - predicted)^2) /
      sum((observed - mean(observed))^2)
  }
  
  for (i in seq_len(n_repeats)) {
    
    iteration_seed <- seed + i - 1
    set.seed(iteration_seed)
    
    train_index <- sample.int(
      n = n,
      size = sample_size,
      replace = FALSE
    )
    
    train_data <- data[train_index, , drop = FALSE]
    test_data <- data[-train_index, , drop = FALSE]
    
    x_train <- as.matrix(
      train_data[, predictors, drop = FALSE]
    )
    
    x_test <- as.matrix(
      test_data[, predictors, drop = FALSE]
    )
    
    dtrain <- xgb.DMatrix(
      data = x_train,
      label = train_data[[response]]
    )
    
    dtest <- xgb.DMatrix(
      data = x_test,
      label = test_data[[response]]
    )
    
    iteration_params <- params
    iteration_params$seed <- iteration_seed
    
    cv_result <- xgb.cv(
      params = iteration_params,
      data = dtrain,
      nrounds = nrounds,
      nfold = nfold,
      early_stopping_rounds = early_stopping_rounds,
      verbose = 0
    )
    
    best_iteration <- cv_result$best_iteration
    
    model <- xgb.train(
      params = iteration_params,
      data = dtrain,
      nrounds = best_iteration,
      verbose = 0
    )
    
    train_prediction <- predict(model, dtrain)
    test_prediction <- predict(model, dtest)
    
    metrics$Best_iteration[i] <- best_iteration
    
    metrics$Train_R2[i] <- calculate_r2(
      train_data[[response]],
      train_prediction
    )
    
    metrics$Test_R2[i] <- calculate_r2(
      test_data[[response]],
      test_prediction
    )
    
    shap_list[[i]] <- shap.values(
      xgb_model = model,
      X_train = x_train
    )$mean_shap_score
  }
  
  shap_matrix <- do.call(cbind, shap_list)
  
  colnames(shap_matrix) <- paste0(
    "Repeat_",
    seq_len(n_repeats)
  )
  
  shap_summary <- data.frame(
    Variable = rownames(shap_matrix),
    Mean_SHAP = rowMeans(shap_matrix),
    SE_SHAP = apply(shap_matrix, 1, sd) /
      sqrt(n_repeats)
  )
  
  shap_summary <- shap_summary[
    order(shap_summary$Mean_SHAP, decreasing = TRUE),
  ]
  
  list(
    metrics = metrics,
    shap_summary = shap_summary,
    shap_matrix = shap_matrix,
    parameters = params,
    predictors = predictors
  )
}

results <- repeat_xgb(
  data = df,
  response = response,
  predictors = predictors,
  params = xgb_params,
  n_repeats = 100,
  nfold = 3,
  sample_size = 50000,
  nrounds = 1000,
  early_stopping_rounds = 30,
  seed = 123
)

# 4. Quantile regression

library(quantreg)
data_raw <- df

variables <- c(
  "Res",
  "Base_Res",
  "AGB",
  "ET",
  "VPD",
  "AI",
  "Clay",
  "Soil_Depth"
)

stopifnot(all(variables %in% names(data_raw)))
data_qr <- data_raw[, variables]
data_qr <- data_qr[complete.cases(data_qr), ]

stopifnot(
  all(vapply(data_qr, is.numeric, logical(1)))
)

winsorize <- function(x, proportion = 0.01) {
  
  limits <- quantile(
    x,
    probs = c(proportion, 1 - proportion),
    na.rm = TRUE,
    type = 7
  )
  
  x[x < limits[1]] <- limits[1]
  x[x > limits[2]] <- limits[2]
  
  x
}

data_qr[] <- lapply(
  data_qr,
  winsorize,
  proportion = 0.01
)

data_qr <- as.data.frame(
  scale(data_qr)
)

quantiles <- seq(
  from = 0.1,
  to = 0.9,
  by = 0.1
)

model_formula <- Res ~
  Base_Res +
  AGB +
  ET + I(ET^2) +
  VPD + I(VPD^2) +
  AI + I(AI^2) +
  Clay +
  Soil_Depth +
  Base_Res:AGB +
  Base_Res:ET

models <- lapply(
  quantiles,
  function(tau_value) {
    rq(
      formula = model_formula,
      tau = tau_value,
      data = data_qr
    )
  }
)


set.seed(123)

coef_list <- lapply(
  seq_along(taus),
  function(i) {
    
    model_summary <- summary(
      rq_models[[i]],
      se = "boot",
      R = 50
    )
    
    coefficients <- as.data.frame(
      model_summary$coefficients
    )
    
    data.frame(
      Quantile = taus[i],
      Variable = rownames(coefficients),
      Estimate = coefficients[, 1],
      Standard_error = coefficients[, 2],
      Test_statistic = coefficients[, 3],
      P_value = coefficients[, 4],
      CI_lower = coefficients[, 1] -
        1.96 * coefficients[, 2],
      CI_upper = coefficients[, 1] +
        1.96 * coefficients[, 2],
      row.names = NULL
    )
  }
)

coefficient_results <- do.call(
  rbind,
  coef_list
)

write.csv(
  coefficient_results,
  "quantile_regression_coefficients.csv",
  row.names = FALSE
)

