library(rms)
library(survival)
library(survivalROC)
library(Hmisc)
library(boot)

data <- read.delim(file.choose(), row.names = 1)

os <- as.numeric(data$OS)
censor <- as.numeric(data$Censor)

S <- Surv(os, censor)
attach(data)

coxm <- cph(Surv(OS, Censor) ~ NK.Riskscore + Path.Riskscore + age + gender + 
              Mx + Nx + Tx + Stage, 
            x = TRUE, y = TRUE, data = data, surv = TRUE)

scoxm <- step(coxm)

factor_vars <- attr(terms(scoxm), "term.labels")
if (length(factor_vars) == 0) {
  factor_vars <- c("NK.Riskscore", "Path.Riskscore", "Tx", "age", "gender", "Stage")
}

factor <- paste(factor_vars, collapse = "+")

ddist <- datadist(data[, factor_vars])
options(datadist = "ddist")

f <- cph(as.formula(paste("Surv(OS, Censor) ~", factor)), 
         x = TRUE, y = TRUE, data = data, surv = TRUE)

surv <- Survival(f)
nom <- nomogram(f, 
                fun = list(function(x) surv(360, x), 
                           function(x) surv(720, x), 
                           function(x) surv(1080, x),
                           function(x) surv(1800, x),
                           function(x) surv(3600, x)), 
                lp = FALSE,
                funlabel = c("1-year survival", "2-year survival", 
                             "3-year survival", "5-year survival", "10-year survival"),
                fun.at = seq(0.90, 0.1, by = -0.1))

plot(nom, cex.axis = 1.2, xfrac = 0.2, lwd = 2, 
     cex.var = 1.2, cex = 1.2, font.axis = 2)

boot_cindex <- function(data, indices, formula) {
  d <- data[indices, ]
  if (sum(d$Censor) < 2) return(0.5)
  fit <- try(cph(formula, data = d, x = TRUE, y = TRUE, surv = TRUE), silent = TRUE)
  if (inherits(fit, "try-error")) return(0.5)
  conc <- concordance(fit)
  return(ifelse(is.na(conc$concordance), 0.5, conc$concordance))
}

c_index_raw <- 1 - as.numeric(rcorrcens(S ~ predict(f)))
c_index_raw <- ifelse(length(c_index_raw) > 1, mean(c_index_raw), c_index_raw)
c_index_raw <- ifelse(is.na(c_index_raw) | is.nan(c_index_raw), 
                      survConcordance(S ~ predict(f))$concordance, c_index_raw)

set.seed(123)
boot_formula <- as.formula(paste("Surv(OS, Censor) ~", factor))
boot_result <- boot(data = data, statistic = boot_cindex, R = 1000, formula = boot_formula)

boot_t <- boot_result$t[boot_result$t >= 0.5 & boot_result$t <= 1]
c_index <- mean(boot_t)
c_se <- sd(boot_t)
c_ci <- try(boot.ci(boot_result, type = "bca")$bca[4:5], silent = TRUE)
if (inherits(c_ci, "try-error") | length(c_ci) < 2) {
  c_ci <- quantile(boot_t, c(0.025, 0.975))
}
c_ci_lower <- max(0.5, c_ci[1])
c_ci_upper <- min(1, c_ci[2])

single_vars <- c("NK.Riskscore", "Path.Riskscore", "Tx", "age", "gender", "Stage")
cindex_list <- list()

for (var in single_vars) {
  ddist_single <- datadist(data[, var])
  options(datadist = "ddist_single")
  f_single <- cph(as.formula(paste("Surv(OS, Censor) ~", var)), 
                  data = data, x = TRUE, y = TRUE, surv = TRUE)
  
  c_single_raw <- 1 - as.numeric(rcorrcens(S ~ predict(f_single)))
  c_single_raw <- ifelse(length(c_single_raw) > 1, mean(c_single_raw), c_single_raw)
  c_single_raw <- ifelse(is.na(c_single_raw) | is.nan(c_single_raw), 
                         survConcordance(S ~ predict(f_single))$concordance, c_single_raw)
  
  boot_single <- boot(data = data, statistic = boot_cindex, R = 1000, 
                      formula = as.formula(paste("Surv(OS, Censor) ~", var)))
  boot_single_t <- boot_single$t[boot_single$t >= 0.5 & boot_single$t <= 1]
  c_single <- mean(boot_single_t)
  c_single_se <- sd(boot_single_t)
  c_single_ci <- try(boot.ci(boot_single, type = "bca")$bca[4:5], silent = TRUE)
  if (inherits(c_single_ci, "try-error") | length(c_single_ci) < 2) {
    c_single_ci <- quantile(boot_single_t, c(0.025, 0.975))
  }
  c_single_ci_lower <- max(0.5, c_single_ci[1])
  c_single_ci_upper <- min(1, c_single_ci[2])
  
  cindex_list[[var]] <- c(c_single, c_single_se, c_single_ci_lower, c_single_ci_upper)
}

cindex_table <- data.frame(
  模型类型 = c("多模态整合模型", "NK.Riskscore单因素", "Path.Riskscore单因素", "Tx单因素", 
           "age单因素", "gender单因素", "Stage单因素"),
  C指数    = round(c(c_index, 
                   cindex_list[["NK.Riskscore"]][1], 
                   cindex_list[["Path.Riskscore"]][1], 
                   cindex_list[["Tx"]][1], 
                   cindex_list[["age"]][1], 
                   cindex_list[["gender"]][1], 
                   cindex_list[["Stage"]][1]), 4),
  标准误   = round(c(c_se, 
                  cindex_list[["NK.Riskscore"]][2], 
                  cindex_list[["Path.Riskscore"]][2], 
                  cindex_list[["Tx"]][2], 
                  cindex_list[["age"]][2], 
                  cindex_list[["gender"]][2], 
                  cindex_list[["Stage"]][2]), 4),
  CI下限   = round(c(c_ci_lower, 
                   cindex_list[["NK.Riskscore"]][3], 
                   cindex_list[["Path.Riskscore"]][3], 
                   cindex_list[["Tx"]][3], 
                   cindex_list[["age"]][3], 
                   cindex_list[["gender"]][3], 
                   cindex_list[["Stage"]][3]), 4),
  CI上限   = round(c(c_ci_upper, 
                   cindex_list[["NK.Riskscore"]][4], 
                   cindex_list[["Path.Riskscore"]][4], 
                   cindex_list[["Tx"]][4], 
                   cindex_list[["age"]][4], 
                   cindex_list[["gender"]][4], 
                   cindex_list[["Stage"]][4]), 4)
)

extract_calibration_data <- function(cal_obj, time_label) {
  cal_mat <- as.matrix(cal_obj)
  cal_df <- data.frame(
    Time     = time_label,
    Predicted = cal_mat[, "mean.predicted"],  
    Observed  = cal_mat[, "KM"],               
    Lower_CI  = cal_mat[, "KM.corrected"],     
    Upper_CI  = cal_mat[, "std.err"]           
  )
  cal_df <- cal_df[!is.na(cal_df$Observed) & is.finite(cal_df$Observed), ]
  return(cal_df)
}

coxm_1y  <- cph(as.formula(paste("Surv(OS, Censor) ~", factor)), x=TRUE, y=TRUE, data=data, surv=TRUE, time.inc=360)
cal_1y   <- calibrate(coxm_1y,  cmethod="KM", method="boot", u=360,  m=40, B=100)
cal1_data <- extract_calibration_data(cal_1y, "1_year")

coxm_2y  <- cph(as.formula(paste("Surv(OS, Censor) ~", factor)), x=TRUE, y=TRUE, data=data, surv=TRUE, time.inc=720)
cal_2y   <- calibrate(coxm_2y,  cmethod="KM", method="boot", u=720,  m=40, B=100)
cal2_data <- extract_calibration_data(cal_2y, "2_year")

coxm_3y  <- cph(as.formula(paste("Surv(OS, Censor) ~", factor)), x=TRUE, y=TRUE, data=data, surv=TRUE, time.inc=1080)
cal_3y   <- calibrate(coxm_3y,  cmethod="KM", method="boot", u=1080, m=40, B=100)
cal3_data <- extract_calibration_data(cal_3y, "3_year")

coxm_5y  <- cph(as.formula(paste("Surv(OS, Censor) ~", factor)), x=TRUE, y=TRUE, data=data, surv=TRUE, time.inc=1800)
cal_5y   <- calibrate(coxm_5y,  cmethod="KM", method="boot", u=1800, m=40, B=100)
cal5_data <- extract_calibration_data(cal_5y, "5_year")

coxm_10y <- cph(as.formula(paste("Surv(OS, Censor) ~", factor)), x=TRUE, y=TRUE, data=data, surv=TRUE, time.inc=3600)
cal_10y  <- calibrate(coxm_10y, cmethod="KM", method="boot", u=3600, m=40, B=100)
cal10_data <- extract_calibration_data(cal_10y, "10_year")

all_cal_data <- rbind(cal1_data, cal2_data, cal3_data, cal5_data, cal10_data)

mae_1y  <- mean(abs(cal1_data$Predicted  - cal1_data$Observed),  na.rm = TRUE)
mae_2y  <- mean(abs(cal2_data$Predicted  - cal2_data$Observed),  na.rm = TRUE)
mae_3y  <- mean(abs(cal3_data$Predicted  - cal3_data$Observed),  na.rm = TRUE)
mae_5y  <- mean(abs(cal5_data$Predicted  - cal5_data$Observed),  na.rm = TRUE)
mae_10y <- mean(abs(cal10_data$Predicted - cal10_data$Observed), na.rm = TRUE)

score_1y  <- 1 - mae_1y
score_2y  <- 1 - mae_2y
score_3y  <- 1 - mae_3y
score_5y  <- 1 - mae_5y
score_10y <- 1 - mae_10y

model_compare_table <- data.frame(
  时间点           = c("1年", "2年", "3年", "5年", "10年"),
  校准误差值_MAE   = round(c(mae_1y, mae_2y, mae_3y, mae_5y, mae_10y), 4),
  校准得分值_1减MAE = round(c(score_1y, score_2y, score_3y, score_5y, score_10y), 4)
)

pred_value <- predict(f)

time_points <- c(360, 720, 1080, 1800, 3600)
time_labels <- c("1年", "2年", "3年", "5年", "10年")
auc_values <- numeric(length(time_points))
names(auc_values) <- time_labels

for (i in seq_along(time_points)) {
  roc_result <- survivalROC(
    Stime       = os, 
    status      = censor, 
    marker      = pred_value,
    predict.time = time_points[i], 
    method      = "KM"
  )
  auc_values[i] <- roc_result$AUC
}

auc_table <- data.frame(
  时间点       = time_labels,
  时间依赖AUC = round(auc_values, 4)
)

model_full_table <- cbind(model_compare_table, 时间依赖AUC = round(auc_values, 4))

roc_1y <- survivalROC(Stime = os, status = censor, marker = pred_value, predict.time = 360,  method = "KM")
roc_2y <- survivalROC(Stime = os, status = censor, marker = pred_value, predict.time = 720,  method = "KM")
roc_3y <- survivalROC(Stime = os, status = censor, marker = pred_value, predict.time = 1080, method = "KM")

plot(roc_1y$FP, roc_1y$TP, type = "l", col = "red", lwd = 3,
     main = "Multimodal COAD Prognostic Model ROC Curve",
     xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1))
lines(roc_2y$FP, roc_2y$TP, col = "blue", lwd = 3, lty = 2)
lines(roc_3y$FP, roc_3y$TP, col = "forestgreen", lwd = 3, lty = 3)
abline(a = 0, b = 1, lty = 2, col = "gray50", lwd = 2)
text(0.7, 0.2, paste("1-year AUC =", round(roc_1y$AUC, 3)), col = "red",    cex = 1.1, font = 2)
text(0.7, 0.1, paste("2-year AUC =", round(roc_2y$AUC, 3)), col = "blue",   cex = 1.1, font = 2)
text(0.7, 0.0, paste("3-year AUC =", round(roc_3y$AUC, 3)), col = "forestgreen", cex = 1.1, font = 2)