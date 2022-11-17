library(dplyr)
library(Seurat)
library(patchwork)
library(Rmisc)

pdf(file="/deac/csc/khuriGrp/zhaok220/data_processing_scRNA/output/plot_byR.pdf")

# Load the PBMC dataset
# 32738 x 2700 sparse Matrix (genes x cells)
# The values in this matrix represent the number of molecules for each feature 
# (i.e. gene; row) that are detected in each cell (column).
pbmc.data <- Read10X(data.dir = "/deac/csc/khuriGrp/zhaok220/data_processing_scRNA/data/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)


# Quality Control
# Filter cells that have unique feature counts over 2,500 or less than 200, and have >5% mitochondrial counts
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") # Low-quality / dying cells often exhibit extensive mitochondrial contamination
plot1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
plot2 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
multiplot(plot1,plot2, cols = 1)

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features(feature selection)
# Highly variable features help highlight biological signal in single-cell datasets.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10) # Identify the 10 most highly variable genes
plot3 <- VariableFeaturePlot(pbmc)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
multiplot(plot3,plot4, cols = 1)

# Scaling
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# PCA reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimHeatmap(pbmc, dims = 1:10, cells = 500, balanced = TRUE) #heatmap of 15 principal components ordered by their PCA scores.

# Cluster cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Non-linear dimensional reduction
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# Identify biomarkers
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
