
library(Matrix)
library(foreach)
library(glmnet)
library(survival)


set.seed(39)


x <- read.delim(file.choose(), row.names = 1) # 读取表达数据
y <- read.delim(file.choose(), row.names = 1) # 读取生存数据


dat <- data.matrix(x)
lable <- Surv(y$OS, y$CENSOR)

cvob1 <- cv.glmnet(dat, lable, family = "cox", alpha = 1)
plot(cvob1)

cvob2 <- glmnet(dat, lable, family = "cox", alpha = 1)
plot(cvob2, xvar = "lambda", ylim = c(-35000000, 35000000), label = TRUE)

coef_min <- coef(cvob1, s = "lambda.min")
coef_1se <- coef(cvob1, s = "lambda.1se")


title("Cox Family", line = 2.5)



result_1se <- data.matrix(coef_1se)
write.table(result_1se, "Lasso result-1se 425 39.txt", quote = FALSE, sep = "\t")

result_min <- data.matrix(coef_min)
write.table(result_min, "Lasso result-min 425 39.txt", quote = FALSE, sep = "\t")


