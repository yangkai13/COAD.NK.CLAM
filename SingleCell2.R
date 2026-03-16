library(dplyr)
library(Seurat)
library(patchwork)
library(celldex)
library(SingleR)
library(BiocParallel)
library(ggplot2)
library(scRNAtoolVis)
library(EnhancedVolcano)
library(plyr)

folder_list <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
seurat_list <- lapply(folder_list, function(folder) {
  data <- Read10X(data.dir = folder)
  CreateSeuratObject(counts = data, project = basename(folder))
})

seurat_merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

seurat_merged[["percent.mt"]] <- PercentageFeatureSet(seurat_merged, pattern = "^MT-")

morandi_colors <- c("#4682B4", "#E59866", "#27AE60", "#C0392B", "#8E44AD", "#3498DB", "#F1C40F")

VlnPlot(seurat_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01, cols = morandi_colors)

FeatureScatter(seurat_merged, "nCount_RNA", "percent.mt", cols = morandi_colors) + 
  FeatureScatter(seurat_merged, "nCount_RNA", "nFeature_RNA", cols = morandi_colors)

seurat_filtered <- subset(seurat_merged, subset = nFeature_RNA > 50 & percent.mt < 5)

seurat <- NormalizeData(seurat_filtered)
seurat <- FindVariableFeatures(seurat, nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)

ElbowPlot(seurat)
seurat_harmony <- IntegrateLayers(object = seurat, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)

seurat_harmony <- FindNeighbors(seurat_harmony, reduction = "harmony", dims = 1:20)
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.5)
seurat_harmony <- RunUMAP(seurat_harmony, reduction = "harmony", dims = 1:20, reduction.name = "umap")

DimPlot(seurat_harmony, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  DimPlot(seurat_harmony, reduction = "umap", group.by = "orig.ident", cols = morandi_colors)

marker_genes <- c("CD3D","CD3E","TRAC","KLRD1","GNLY","NKG7","CD19","CD79A","MS4A1",
                  "CD68","CD163","CLEC9A","XCR1","KIT","TPSAB1","TPSB2",
                  "COL1A1","COL3A1","VWF","PECAM1","EPCAM","KRT18","KRT19",
                  "FGFBP2","FCGR3A","CX3CR1","PRF1")

VlnPlot(seurat_harmony, features = marker_genes, ncol = 4)
FeaturePlot(seurat_harmony, features = marker_genes, reduction = "umap")

seurat_harmony$celltype <- plyr::mapvalues(
  seurat_harmony$seurat_clusters,
  from = 0:14,
  to = c("T_cell","T_cell","T_cell","B_cell","B_cell","B_cell",
         "Myeloid","NK_cell","Epithelial","B_cell","B_cell",
         "Fibroblast","Myeloid","Mast","T_cell")
)

DimPlot(seurat_harmony, reduction = "umap", group.by = "celltype", label = TRUE)

Idents(seurat_harmony) <- "celltype"
deg_all <- FindAllMarkers(seurat_harmony, min.pct = 0.25, logfc.threshold = 0.25)

jjVolcano(diffData = deg_all, tile.col = morandi_colors, myMarkers = c("CD3D","EPCAM","CD79A","COL1A1","GNLY","NKG7"))

prop_data <- as.data.frame(table(seurat_harmony$orig.ident, seurat_harmony$celltype))
colnames(prop_data) <- c("sample", "celltype", "count")

ggplot(prop_data, aes(x = sample, weight = count, fill = celltype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = morandi_colors) +
  theme_minimal() +
  RotatedAxis() +
  labs(y = "Proportion")

nk <- subset(seurat_harmony, celltype == "NK_cell")

nk <- FindVariableFeatures(nk)
nk <- ScaleData(nk)
nk <- RunPCA(nk)
nk <- FindNeighbors(nk, dims = 1:20)
nk <- FindClusters(nk, resolution = 0.5)
nk <- RunUMAP(nk, dims = 1:20)

DimPlot(nk, label = TRUE)

nk_markers_key <- c("NCAM1","FCGR3A","CX3CR1","IL7R","CCL5","IL32","CREM","PRF1","GNLY","NKG7")
VlnPlot(nk, features = nk_markers_key)
FeaturePlot(nk, features = nk_markers_key)

nk$nk_subtype <- plyr::mapvalues(
  nk$seurat_clusters,
  from = 0:6,
  to = c("NK2","NK2","NK1","NK2","NK2","NK1","NK2")
)

DimPlot(nk, group.by = "nk_subtype", label = TRUE)

nk_diff <- FindMarkers(
  nk,
  ident.1 = c("GSM4904237_Patient1_normal", "GSM4904240_Patient2_normal", "GSM4904246_Patient3_normal"),
  ident.2 = c("GSM4904235_Patient1_adenoma", "GSM4904238_Patient2_adenoma", 
              "GSM4904242_Patient3_adenoma", "GSM4904243_Patient3_adenoma_2"),
  group.by = "orig.ident",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

EnhancedVolcano(
  nk_diff,
  lab = rownames(nk_diff),
  x = "avg_log2FC",
  y = "p_val_adj",
  pCutoff = 0.05,
  FCcutoff = 0.5
)