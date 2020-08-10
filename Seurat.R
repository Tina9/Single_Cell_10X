library(Seurat)

library(dplyr)
library(patchwork)

setwd("/Users/xuzhang/Dropbox/Share with Xu/Latest Version/Seurat/MERGE_1_2/")
pbmc.data <- Read10X(data.dir = "/Users/xuzhang/Dropbox/Share with Xu/Latest Version/Seurat/MERGE_1_2/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "sample2",
                           min.cells = 3, min.features = 200)


#### Standard pre-processing workflow
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

## Visualizing QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


## filter conditions: 1.unique feature counts over 2500 or less than 200 2. >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 5)

#### Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1000) ### the same as: pbmc <- NormalizeData(pbmc)

#### Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 20)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
plot2

## Scaling the data
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc, features = all.gens)

# #### Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# print(pbmc[["pca"]], dims = 1:20, nfeatures = 20)
# VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# DimPlot(pbmc, reduction = "pca")
# DimHeatmap(pbmc, dims = 1, cells = 500, balanced = T)
# DimHeatmap(pbmc, dims = 1:20, cells = 500, balanced = T)

#### Determine the 'dimensionality' of the dataset
## Using resampling test
# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# JackStrawPlot(pbmc, dims = 1:15)
# 
# ## Using ElbowPlot function
# ElbowPlot(pbmc, ndims = 50)

#### Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#### Run non-linear dimensional reduction (UMAP/tSNE)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = T)
# Save the object at this point
saveRDS(pbmc, file = "../pbmc_tutoria.rds")

#### Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = T)
VlnPlot(pbmc, features = c("RPS27", "S100A9"))

# Plotting raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = T)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A",
                               "FCGR3A", "LYZ", "PPBP", "CD8A"))

# Generating an expression heatmap for given cells and features
top10 <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#### Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono",
                     "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
saveRDS(pbmc, file = "../pbmc3k_final.rds")
