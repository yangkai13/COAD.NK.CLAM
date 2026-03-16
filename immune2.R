# 脚本功能：基于TCGA COAD数据集（TPM格式）进行基因表达数据预处理，运行CIBERSORT算法推断细胞组分，并生成可视化结果
# 依赖的R包：
#   - e1071：支持向量机算法，用于CIBERSORT核心计算
#   - preprocessCore：数据标准化处理
#   - limma：生物信息学数据分析，包含avereps函数
#   - ggplot2、tidyverse、reshape2：数据可视化
# 运行前需确保安装以下包：
# install.packages(c("e1071", "tidyverse", "ggplot2"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("preprocessCore", "limma"))

# 加载R包
library(e1071)
library(preprocessCore)
library(limma)
library(ggplot2)
library(tidyverse)
library(reshape2)

# 设置工作目录和输入文件
setwd("/Users/q")  # 设置工作目录，请替换为实际路径
expFile <- "TCGA.COAD.TPM.gm.txt"  # 输入文件：TCGA COAD基因表达数据（TPM格式）

# 1. 数据预处理
# 读取基因表达数据文件，文件以制表符分隔，包含表头
rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)  # 转换为矩阵
rownames(rt) <- rt[, 1]  # 第一列设为行名（基因名称）
exp <- rt[, 2:ncol(rt)]  # 提取表达数据（去除第一列）

# 设置矩阵维度名称
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)

# 平均重复基因的表达值（使用limma包的avereps函数）
data <- avereps(data)

# 过滤掉行均值为0的基因，保留有表达的基因
data <- data[rowMeans(data) > 0, ]

# 2. 数据保存
# 直接使用处理后的数据（TPM格式，无需voom转换）
out <- data
out <- rbind(ID = colnames(out), out)  # 添加样本ID作为第一行
# 保存处理后的数据到文件
write.table(out, file = "uniq.symbol.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# 3. 运行CIBERSORT分析
# 加载CIBERSORT脚本（需确保CIBERSORT.R文件在工作目录中）
source("CIBERSORT.R")

# 运行CIBERSORT算法
# 参数说明：
#   - ref.txt：参考基因表达文件（细胞特异性签名矩阵）
#   - uniq.symbol.txt：处理后的混合样本基因表达数据
#   - perm=10：置换检验次数（设为10，实际分析建议>=100）
#   - QN=FALSE：不进行分位数标准化（TPM数据推荐设置为FALSE）
results <- CIBERSORT("ref.txt", "uniq.symbol.txt", perm = 10, QN = FALSE)

# 保存CIBERSORT结果到文件
write.table(results, file = "res_cibersort.txt", sep = "\t", col.names = TRUE, 
            row.names = TRUE, quote = FALSE)

# 4. 可视化结果
# 定义莫兰迪色系（对比度强，适合区分细胞类型）
colour <- c(
  "#B7C3B6", "#D8A7B1", "#8B9A8F", "#E8B5AC", "#6D8299", "#C4B7CB",
  "#A8C4B8", "#F4C4B1", "#7B8BA8", "#D9C2D6", "#9BB8A5", "#E8A9A2",
  "#6A8296", "#C8B6C8", "#A2C1B0", "#F2B5A7", "#7C8EA5", "#D4B8C4",
  "#8DB5A0", "#E8B0A5", "#6B8799", "#CBB4C9", "#A5C0B3"
)  # 23种颜色，覆盖22种细胞类型+1备用

# 定义自定义ggplot主题
my_theme <- function() {
  theme(
    panel.grid = element_blank(),  # 去除网格线
    panel.border = element_blank(),  # 去除面板边框
    legend.position = "right",  # 图例位置
    legend.text = element_text(size = 8),  # 图例文本大小
    legend.title = element_text(size = 8),  # 图例标题大小
    axis.line = element_line(size = 1),  # 坐标轴线粗细
    text = element_text(family = "Times"),  # 文本字体
    axis.text.y = element_text(size = 8, face = "bold", color = "black"),  # y轴标签样式
    axis.text.x = element_text(size = 8, face = "bold", color = "black", 
                               angle = 90, hjust = 1),  # x轴标签样式，旋转90度
    axis.title = element_text(size = 10, face = "bold"),  # 坐标轴标题样式
    plot.title = element_text(hjust = 0.5, size = 10)  # 图形标题居中
  )
}

# 准备绘图数据：提取前22列（细胞类型比例），转换为长格式
p1 <- results[, 1:22] %>% 
  reshape2::melt() %>%  # 转换为长格式
  ggplot(aes(x = Var1, y = value, fill = Var2)) +  # 绘制堆叠柱状图
  geom_bar(stat = "identity") +  # 堆叠柱状图
  coord_flip() +  # 翻转坐标轴
  scale_fill_manual(values = colour) +  # 使用莫兰迪色系
  theme_bw() +  # 使用简洁主题
  theme(
    panel.border = element_blank(),  # 去除面板边框
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank(),  # 去除次网格线
    axis.line = element_line(colour = "black")  # 设置坐标轴线颜色
  ) + 
  my_theme()  # 应用自定义主题

# 保存图形到PDF文件
pdf("cibersort.pdf")
print(p1)
dev.off()
