# -----------------------------------
# Single-Cell RNA-seq Analysis (PBMC)
# Full Pipeline with Seurat in R
# -----------------------------------

### 1. Install and Load Packages ###
install.packages("Seurat")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("patchwork")  # For combining plots

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

### 2. Download and Load Data ###
# Note: Manually download from:
# https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
# Extract to folder "pbmc3k_filtered_gene_bc_matrices"

# Load the data
pbmc_data <- Read10X(data.dir = "C:/Users/adwik/Downloads/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data, project = "pbmc3k")

### 3. Quality Control (QC) ###
# Calculate mitochondrial gene percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
qc_plots <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    ncol = 3, pt.size = 0.1)
print(qc_plots)

# Filter cells:
# - Keep cells with >200 genes (nFeature_RNA)
# - Keep cells with <2500 genes (remove doublets)
# - Keep cells with <5% mitochondrial genes
pbmc <- subset(pbmc, subset = 
                 nFeature_RNA > 200 & 
                 nFeature_RNA < 2500 & 
                 percent.mt < 5)

### 4. Normalization ###
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

### 5. Feature Selection ###
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# Plot variable features
var_feat_plot <- VariableFeaturePlot(pbmc)
var_feat_plot <- LabelPoints(plot = var_feat_plot, points = top10, repel = TRUE)
print(var_feat_plot)

### 6. Scaling and PCA ###
all_genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all_genes)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

# Elbow plot to determine significant PCs
ElbowPlot(pbmc)

### 7. Clustering ###
pbmc <- FindNeighbors(pbmc, dims = 1:10)  # Using first 10 PCs
pbmc <- FindClusters(pbmc, resolution = 0.5)

### 8. UMAP Visualization ###
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)

### 9. Marker Identification ###
# Find markers for every cluster compared to all remaining cells
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top 2 markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

print(top_markers)

### 10. Cell Type Annotation ###
# Based on known markers:
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", 
                     "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Plot annotated UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

### 11. Save Results ###
saveRDS(pbmc, file = "pbmc_processed.rds")

