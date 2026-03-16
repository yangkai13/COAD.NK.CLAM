
# 加载必要的R包
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(plyr)

# 设置工作目录（请根据实际情况修改）
setwd("")

# 1.1 导入数据
folder_list <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
seurat_list <- list()

for (folder in folder_list) {
  data <- Read10X(data.dir = folder)
  seurat_obj <- CreateSeuratObject(counts = data, project = basename(folder))
  seurat_list[[basename(folder)]] <- seurat_obj
}

# 合并所有样本
seurat_merge_COAD_normal <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)], 
                                  add.cell.ids = names(seurat_list))
saveRDS(seurat_merge_COAD_normal, file = "1.最初融合seurat_merge_COAD_normal.rds")

seurat_merge_COAD_normal <- readRDS("1.最初融合seurat_merge_COAD_normal.rds")

# 1.2 质量控制（QC） - 计算线粒体比例
seurat_merge_COAD_normal[["percent.mt"]] <- PercentageFeatureSet(seurat_merge_COAD_normal, pattern = "^MT-")

# 定义颜色方案
morandi_colors <- c("#4682B4", "#E59866", "#27AE60", "#C0392B", "#8E44AD", "#3498DB", "#F1C40F")

# 查看原始细胞总数
ncells_raw <- ncol(seurat_merge_COAD_normal)
cat("合并后的总细胞数（原始）：", ncells_raw, "\n")  # 应为 30896

# 保存过滤前QC图（补充图素材）
p_vln_raw <- VlnPlot(seurat_merge_COAD_normal, 
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                     ncol = 3, pt.size = 0, cols = morandi_colors) + 
  ggtitle("Pre-filtering QC metrics")
ggsave("QC_pre_filter_violin.png", p_vln_raw, width = 10, height = 4)

p_scatter_raw1 <- FeatureScatter(seurat_merge_COAD_normal, "nCount_RNA", "percent.mt", cols = morandi_colors)
p_scatter_raw2 <- FeatureScatter(seurat_merge_COAD_normal, "nCount_RNA", "nFeature_RNA", cols = morandi_colors)
p_scatter_raw <- p_scatter_raw1 + p_scatter_raw2 + plot_annotation(title = "Pre-filtering Scatter Plots")
ggsave("QC_pre_filter_scatter.png", p_scatter_raw, width = 12, height = 5)

# ────────────────────────────────────────────────
# 分支1：宽松阈值过滤（主流程）
# ────────────────────────────────────────────────

seurat_permissive <- subset(seurat_merge_COAD_normal, subset = nFeature_RNA > 50 & percent.mt < 5)

ncells_permissive <- ncol(seurat_permissive)
cat("宽松阈值过滤后细胞数：", ncells_permissive, "\n")  # 应为 12042

# 计算宽松阈值QC统计（防错写法）
meta_perm <- seurat_permissive@meta.data
qc_stats_permissive <- data.frame(
  median_nFeature     = median(meta_perm$nFeature_RNA, na.rm = TRUE),
  median_nCount       = median(meta_perm$nCount_RNA,   na.rm = TRUE),
  median_percent_mt   = median(meta_perm$percent.mt,   na.rm = TRUE),
  n_cells             = nrow(meta_perm)
)
cat("宽松阈值QC统计：\n")
print(qc_stats_permissive)

# 保存宽松阈值后QC图
p_vln_perm <- VlnPlot(seurat_permissive, 
                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                      ncol = 3, pt.size = 0, cols = morandi_colors) + 
  ggtitle("Post-filtering (permissive: nFeature >50, mt<5%)")
ggsave("QC_post_permissive_violin.png", p_vln_perm, width = 20, height = 6)

# 主流程标准化 & 高变基因 & PCA & Harmony（精简版）
seurat_permissive <- NormalizeData(seurat_permissive, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_permissive <- FindVariableFeatures(seurat_permissive, selection.method = "vst", nfeatures = 2000)
seurat_permissive <- ScaleData(seurat_permissive)
seurat_permissive <- RunPCA(seurat_permissive, features = VariableFeatures(seurat_permissive))
seurat_permissive <- FindNeighbors(seurat_permissive, dims = 1:20)
seurat_permissive <- FindClusters(seurat_permissive, resolution = 0.5)
seurat_permissive <- RunUMAP(seurat_permissive, dims = 1:20)

# （此处可继续原代码的Harmony批次校正、细胞注释、NK亚群分析等）
# saveRDS(seurat_permissive, "main_permissive_harmony.rds")  # 保存主流程对象

# ────────────────────────────────────────────────
# 分支2：严格阈值过滤（完整敏感性分析）
# ────────────────────────────────────────────────

seurat_strict <- subset(seurat_merge_COAD_normal, subset = nFeature_RNA > 200 & percent.mt < 5)

ncells_strict <- ncol(seurat_strict)
cat("严格阈值过滤后细胞数：", ncells_strict, "\n")

# 计算严格阈值QC统计
meta_strict <- seurat_strict@meta.data
qc_stats_strict <- data.frame(
  median_nFeature     = median(meta_strict$nFeature_RNA, na.rm = TRUE),
  median_nCount       = median(meta_strict$nCount_RNA,   na.rm = TRUE),
  median_percent_mt   = median(meta_strict$percent.mt,   na.rm = TRUE),
  n_cells             = nrow(meta_strict)
)
cat("严格阈值QC统计：\n")
print(qc_stats_strict)

# 保存严格阈值后QC图
p_vln_strict <- VlnPlot(seurat_strict, 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                        ncol = 3, pt.size = 0, cols = morandi_colors) + 
  ggtitle("Post-filtering (strict: nFeature >200, mt<5%)")
ggsave("QC_post_strict_violin.png", p_vln_strict, width = 10, height = 4)

# 严格阈值完整下游流程（与主流程相同步骤）
seurat_strict <- NormalizeData(seurat_strict, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_strict <- FindVariableFeatures(seurat_strict, selection.method = "vst", nfeatures = 2000)
seurat_strict <- ScaleData(seurat_strict)
seurat_strict <- RunPCA(seurat_strict, features = VariableFeatures(seurat_strict))
seurat_strict <- FindNeighbors(seurat_strict, dims = 1:20)
seurat_strict <- FindClusters(seurat_strict, resolution = 0.5)
seurat_strict <- RunUMAP(seurat_strict, dims = 1:20)

# 保存严格阈值对象（用于对比）
saveRDS(seurat_strict, "strict_threshold_full_pipeline.rds")


p_umap_perm <- DimPlot(seurat_permissive, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("Permissive threshold UMAP")
p_umap_strict <- DimPlot(seurat_strict, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("Strict threshold UMAP")

p_compare <- p_umap_perm + p_umap_strict + plot_annotation(title = "Threshold Comparison: UMAP")
ggsave("UMAP_threshold_comparison.png", p_compare, width = 14, height = 6)

cat("分析完成。已生成QC图和阈值对比图，可用于补充材料。\n")
cat("宽松阈值QC统计请见控制台输出；严格阈值同理。\n")