library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
S1.data <- Read10X(data.dir = "~/Downloads/Sample1/")
S1 <- CreateSeuratObject(counts = S1.data, project = "S1", min.cells = 3, min.features = 200)

# Load the PBMC dataset
S2.data <- Read10X(data.dir = "~/Downloads/Sample2/")
S2 <- CreateSeuratObject(counts = S2.data, project = "S2", min.cells = 3, min.features = 200)

pancreas.list <- list()
pancreas.list[["S1"]] <- NormalizeData(S1, verbose = FALSE)
pancreas.list[["S1"]] <- FindVariableFeatures(S1, selection.method = "vst", 
                                           nfeatures = 2000, verbose = FALSE)

pancreas.list[["S2"]] <- NormalizeData(S2, verbose = FALSE)
pancreas.list[["S2"]] <- FindVariableFeatures(S2, selection.method = "vst", 
                                              nfeatures = 2000, verbose = FALSE)

reference.list <- pancreas.list[c("S1", "S2")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)


library(ggplot2)
library(cowplot)
library(patchwork)
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
pancreas.integrated <- RunTSNE(pancreas.integrated, dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(pancreas.integrated, reduction = "pca", group.by = "orig.ident")
p3 <- DimPlot(pancreas.integrated, reduction = "tsne", group.by = "orig.ident")
p1 + p2

ElbowPlot(pancreas.integrated, ndims = 50)
pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:30)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)
p4 <- DimPlot(pancreas.integrated, reduction = "umap")

pancreas.integrated <- RunTSNE(pancreas.integrated, dims = 1:30)
DimPlot(pancreas.integrated, reduction = "tsne")

pancreas.integrated.markers <- FindAllMarkers(pancreas.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas.integrated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

p5 <- FeaturePlot(pancreas.integrated,
            features = c("S100a9", "Cd14", "Fcgr3",
                         "Nlrp3", "C1qc", "Cdc45"))
p4 + p5

top10 <- pancreas.integrated.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(pancreas.integrated, features = top10$gene)
