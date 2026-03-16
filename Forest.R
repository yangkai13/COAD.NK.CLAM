library(randomForest)
library(randomForestSRC)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(survAUC)
library(pROC)

survival_data <- read.delim("OS.txt", sep = "\t", stringsAsFactors = FALSE)
gene_data     <- read.delim("54维.txt", sep = "\t", stringsAsFactors = FALSE)

colnames(survival_data)[colnames(survival_data) == "ID"] <- "sample_id"
colnames(gene_data)[colnames(gene_data) == "ID"]         <- "sample_id"

merged_data <- merge(survival_data, gene_data, by = "sample_id")

X <- merged_data %>% select(-sample_id, -CENSOR, -OS)
y <- Surv(time = merged_data$OS, event = merged_data$CENSOR)

X_scaled <- as.data.frame(scale(X))
data_for_model <- cbind(X_scaled, OS = merged_data$OS, CENSOR = merged_data$CENSOR)

set.seed(42)
train_indices <- sample(1:nrow(data_for_model), 0.7 * nrow(data_for_model))
train_data <- data_for_model[train_indices, ]
test_data  <- data_for_model[-train_indices, ]

rsf_model <- rfsrc(Surv(OS, CENSOR) ~ ., 
                   data       = train_data,
                   ntree      = 500, 
                   mtry       = ncol(X_scaled) / 3, 
                   nodesize   = 5, 
                   splitrule  = "logrank",
                   importance = TRUE,
                   seed       = 42)

rsf_pred_test <- predict(rsf_model, newdata = test_data)
c_index_test  <- rsf_pred_test$err.rate[length(rsf_pred_test$err.rate)]

feature_importance <- data.frame(
  Feature    = names(rsf_model$importance),
  Importance = rsf_model$importance
) %>% arrange(desc(Importance))

rsf_pred_all <- predict(rsf_model, newdata = data_for_model)
risk_scores_all <- rsf_pred_all$predicted

median_threshold <- median(risk_scores_all)

data_for_model$Risk_Group <- ifelse(risk_scores_all > median_threshold, "High Risk", "Low Risk")

surv_obj_all <- Surv(time = data_for_model$OS, event = data_for_model$CENSOR)
km_fit_all   <- survfit(surv_obj_all ~ Risk_Group, data = data_for_model)

ggsurv_all <- ggsurvplot(km_fit_all,
                         data              = data_for_model, 
                         pval              = TRUE, 
                         pval.method       = TRUE,
                         pval.size         = 6,
                         risk.table        = TRUE,
                         risk.table.fontsize = 6,
                         risk.table.height = 0.2,
                         title             = "Kaplan-Meier Survival Curves by Risk Group (All Samples)",
                         xlab              = "Time (days)",
                         ylab              = "Survival Probability",
                         legend.title      = "Risk Group",
                         legend.labs       = c("High Risk", "Low Risk"),
                         palette           = c("#D73027", "#4575B4"),
                         font.main         = c(14, "bold", "black"),
                         font.x            = c(20, "plain", "black"),
                         font.y            = c(20, "plain", "black"),
                         font.tickslab     = c(14, "plain", "black"),
                         font.legend       = c(14, "plain", "black"),
                         conf.int          = TRUE)

arrange_ggsurvplots(list(ggsurv_all), print = FALSE, nrow = 1, ncol = 1)

risk_scores_test <- rsf_pred_test$predicted
test_data$Risk_Group <- ifelse(risk_scores_test > median_threshold, "High Risk", "Low Risk")

surv_obj_test <- Surv(time = test_data$OS, event = test_data$CENSOR)
km_fit_test   <- survfit(surv_obj_test ~ Risk_Group, data = test_data)

ggsurv_test <- ggsurvplot(km_fit_test,
                          data              = test_data,
                          pval              = TRUE,
                          pval.method       = TRUE,
                          pval.size         = 6,
                          risk.table        = TRUE,
                          risk.table.fontsize = 6,
                          risk.table.height = 0.2,
                          title             = "Kaplan-Meier Survival Curves by Risk Group (Test Set)",
                          xlab              = "Time (days)",
                          ylab              = "Survival Probability",
                          legend.title      = "Risk Group",
                          legend.labs       = c("High Risk", "Low Risk"),
                          palette           = c("#D73027", "#4575B4"),
                          font.main         = c(14, "bold", "black"),
                          font.x            = c(20, "plain", "black"),
                          font.y            = c(20, "plain", "black"),
                          font.tickslab     = c(14, "plain", "black"),
                          font.legend       = c(14, "plain", "black"),
                          conf.int          = TRUE)

arrange_ggsurvplots(list(ggsurv_test), print = FALSE, nrow = 1, ncol = 1)

roc_obj_all <- roc(response = data_for_model$CENSOR, predictor = risk_scores_all,
                   levels = c(0, 1), direction = "<")

roc_obj_train <- roc(response = train_data$CENSOR, predictor = rsf_pred_train$predicted,
                     levels = c(0, 1), direction = "<")

roc_obj_test  <- roc(response = test_data$CENSOR,  predictor = risk_scores_test,
                     levels = c(0, 1), direction = "<")

class_data <- cbind(sample_id = merged_data$sample_id, X, group = data_for_model$Risk_Group)

set.seed(123456)
rf <- randomForest(as.factor(group) ~ ., data = class_data %>% select(-sample_id), ntree = 500)

optionTrees <- which.min(rf$err.rate[,1])
rf2 <- randomForest(as.factor(group) ~ ., data = class_data %>% select(-sample_id), ntree = optionTrees)

importance <- importance(rf2)
importance_df <- data.frame(
  ID         = names(importance[, "MeanDecreaseGini"]),
  Importance = importance[, "MeanDecreaseGini"]
) %>% arrange(desc(Importance))

rfGenes <- names(importance[order(importance[, "MeanDecreaseGini"], decreasing = TRUE),])[importance[order(importance[, "MeanDecreaseGini"], decreasing = TRUE),] > 2]

sigExp <- t(X[, rfGenes, drop = FALSE])
sigExpOut <- rbind(ID = colnames(sigExp), sigExp)