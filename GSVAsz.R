# 加载必要的R包
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)

# 文件路径和参数设置
expFile <- "merge.txt"               # 表达数据文件
clusterFile <- "PRGcluster.txt"      # 聚类结果文件
gmtFile <- "c5.go.v7.5.1.symbols.txt"          # 基因集文件
setwd("/Users/q")                    # 设置工作目录

# 自定义颜色设置
heatmap_colors <- colorRampPalette(c("cadetblue", "white", "red4"))(50)  # 热图颜色渐变
ann_colors <- list(
  PRGcluster = c("L" = "#0066FF", "H" = "#FF9900")  # 为L和H聚类定义颜色
)

# 读取表达数据文件并预处理
rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)  # 去重处理，处理重复基因

# GSVA分析
geneSets <- getGmt(gmtFile, geneIdType = SymbolIdentifier())
gsvaResult <- gsva(data, 
                   geneSets, 
                   min.sz = 10, 
                   max.sz = 500, 
                   verbose = TRUE, 
                   parallel.sz = 1)
# 输出GSVA结果
gsvaOut <- rbind(id = colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file = "gsvaOut.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# 读取聚类文件
cluster <- read.table(clusterFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# 数据合并
gsvaResult <- t(gsvaResult)
sameSample <- intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult <- gsvaResult[sameSample, , drop = FALSE]
cluster <- cluster[sameSample, , drop = FALSE]
gsvaCluster <- cbind(gsvaResult, cluster)
Project <- gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
gsvaCluster <- cbind(gsvaCluster, Project)

# 差异分析
adj.P.Val.Filter <- 0.05  # 调整P值阈值
allType <- as.vector(gsvaCluster$PRGcluster)
comp <- combn(levels(factor(allType)), 2)  # 生成所有可能的聚类对比较（L vs H）

for (i in 1:ncol(comp)) {
  # 按聚类分组
  treat <- gsvaCluster[gsvaCluster$PRGcluster == comp[2, i], ]
  con <- gsvaCluster[gsvaCluster$PRGcluster == comp[1, i], ]
  data <- rbind(con, treat)
  
  # 构建通路差异分析
  Type <- as.vector(data$PRGcluster)
  ann <- data[, c(ncol(data), (ncol(data) - 1))]  # 提取Project和PRGcluster列作为注释
  data <- t(data[, -c((ncol(data) - 1), ncol(data))])  # 移除注释列，保留GSVA数据
  design <- model.matrix(~0 + factor(Type))
  colnames(design) <- levels(factor(Type))
  fit <- lmFit(data, design)
  contrast <- paste0(comp[2, i], "-", comp[1, i])  # 构建对比，如"H-L"
  cont.matrix <- makeContrasts(contrast, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  # 输出所有通路差异结果
  allDiff <- topTable(fit2, adjust = 'fdr', number = 200000)
  allDiffOut <- rbind(id = colnames(allDiff), allDiff)
  write.table(allDiffOut, file = paste0(contrast, ".all.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
  
  # 筛选显著差异通路
  diffSig <- allDiff[with(allDiff, (abs(logFC) > 0.1 & adj.P.Val < adj.P.Val.Filter)), ]
  diffSigOut <- rbind(id = colnames(diffSig), diffSig)
  write.table(diffSigOut, file = paste0(contrast, ".diff.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
  
  # 动态生成注释颜色
  current_clusters <- c(comp[1, i], comp[2, i])
  ann_colors_subset <- list(PRGcluster = ann_colors$PRGcluster[current_clusters])
  
  # 绘制差异通路热图
  termNum <- 20  # 热图显示的通路数量
  diffTermName <- as.vector(rownames(diffSig))
  diffLength <- length(diffTermName)
  if (diffLength < termNum) termNum <- diffLength
  hmGene <- diffTermName[1:termNum]
  hmExp <- data[hmGene, ]
  
  pdf(file = paste0(contrast, ".heatmap.pdf"), width = 16, height = 6)
  pheatmap(hmExp, 
           annotation = ann,
           annotation_colors = ann_colors_subset,
           color = heatmap_colors,  # 使用自定义热图颜色
           cluster_cols = FALSE,
           show_colnames = FALSE,
           gaps_col = as.vector(cumsum(table(Type))),  # 按聚类分组分隔
           scale = "row",
           fontsize = 8,
           fontsize_row = 6,
           fontsize_col = 8)
  dev.off()
}

