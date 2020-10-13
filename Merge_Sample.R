library(Seurat)
library(dplyr)
library(patchwork)
library(cowplot)
library(ggplot2)
library(SingleR)
library(reshape2)
library(RColorBrewer)

HP4540PR <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4540_1"
HP4540RE <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/aggr2_3/"

HP4568PR <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4568Pr/"
HP4568RE <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4568Re/"

HP4313PR <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4313Pr/"
HP4313RE <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4313Re/"

HP4540PR.data <- Read10X(data.dir = HP4540PR)
HP4540PR <- CreateSeuratObject(counts = HP4540PR.data, project = "HP4540PR")

HP4540RE.data <- Read10X(data.dir = HP4540RE)
HP4540RE <- CreateSeuratObject(counts = HP4540RE.data, project = "HP4540RE")

HP4568PR.data <- Read10X(data.dir = HP4568PR)
HP4568PR <- CreateSeuratObject(counts = HP4568PR.data, project = "HP4568PR")

HP4568RE.data <- Read10X(data.dir = HP4568RE)
HP4568RE <- CreateSeuratObject(counts = HP4568RE.data, project = "HP4568RE")

HP4313PR.data <- Read10X(data.dir = HP4313PR)
HP4313PR <- CreateSeuratObject(counts = HP4313PR.data, project = "HP4313PR")

HP4313RE.data <- Read10X(data.dir = HP4313RE)
HP4313RE <- CreateSeuratObject(counts = HP4313RE.data, project = "HP4313RE")

pbmc.combined <- merge(HP4540PR,
                       y = c(HP4540RE, HP4568PR, HP4568RE, HP4313PR, HP4313RE),
                       add.cell.ids = c("HP4540PR", "HP4540RE", "HP4568PR", "HP4568RE", "HP4313PR", "HP4313RE"),
                       project = "PARPi")

pbmc.combined[["percent.mt"]] <- PercentageFeatureSet(pbmc.combined, pattern = "^mt-")
#density_plot(pbmc.combined@meta.data)
pbmc.combined <- subset(pbmc.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 5)
pbmc.combined <- NormalizeData(pbmc.combined)

######## Finding Markers #############
pbmc.combined <- FindVariableFeatures(pbmc.combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc.combined)
pbmc.combined <- ScaleData(pbmc.combined)
pbmc.combined <- RunPCA(pbmc.combined, npcs = 30, features = VariableFeatures(object = pbmc.combined))
DimHeatmap(pbmc.combined, dims = 1, cells = 500, balanced = T)
ElbowPlot(pbmc.combined, ndims = 30)

pbmc.combined <- FindNeighbors(pbmc.combined, dims = 1:21)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.5)
pbmc.combined <- RunUMAP(pbmc.combined, dims = 1:21)

p1 <- DimPlot(pbmc.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(pbmc.combined, reduction = "umap", label = T) 
p3 <- DimPlot(pbmc.combined, reduction = "umap", label = T, split.by = "orig.ident") + NoLegend()
top_row <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
p4 <- plot_grid(top_row, p3, labels = c('', 'C'), label_size = 12, ncol = 1)

pbmc.combined <- RunTSNE(pbmc.combined, dims = 1:21)
p1 <- DimPlot(pbmc.combined, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(pbmc.combined, reduction = "tsne", label = T) 
p3 <- DimPlot(pbmc.combined, reduction = "tsne", label = T, split.by = "orig.ident") + NoLegend()
top_row <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
p4 <- plot_grid(top_row, p3, labels = c('', 'C'), label_size = 12, ncol = 1)
plot(p4)

DimPlot(pbmc.combined, reduction = "pca", group.by = "orig.ident")

####### Seeking for clusters
all_markers <- FindAllMarkers(pbmc.combined, only.pos = 0.25, logfc.threshold = 0.25)
final_markers <- all_markers %>% group_by(cluster) %>% as.data.frame()

top10 <- all_markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_logFC)

metadata <- pbmc.combined@meta.data %>% as.data.frame
########## Validation Clusters
mouse_refer <- ImmGenData()
ExpressData <- GetAssayData(pbmc.combined, slot = "data")  
annote_res <- SingleR(test = ExpressData, ref = mouse_refer, labels = mouse_refer$label.main)
metadata$Immlabels <- as.data.frame(annote_res)$labels


# write.table(metadata, file = "~/Desktop/merge_metadata.txt", sep = "\t", quote = F)
# write.table(final_markers[final_markers$p_val < 0.05, ], file = "~/Desktop/marker_genes.txt",
#             sep = "\t", quote = F)

####### Statisticing the change of clusters ##############
component_stat <- function(metadata, clus_num){
  metadata$samples <- sapply(strsplit(rownames(metadata), split = "_"), "[[", 1)
  plot_stat <- metadata %>% 
    group_by(orig.ident) %>%
    dplyr::count(seurat_clusters) %>%
    as.data.frame()
  plot_stat$sample <- substr(plot_stat$orig.ident, 1, 6)
  plot_percent_stat <- plot_stat %>%
    group_by(orig.ident) %>%
    mutate(per = (n/sum(n) * 100))
  
  nb.cols <- clus_num
  mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)
  
  ggplot(data = plot_percent_stat, aes(x = orig.ident, y = per, fill = seurat_clusters)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = mycolors) +
    #facet_grid(~sample) +
    #theme_classic() +
    xlab("Samples") +
    ylab("Cluster Percent") +
    labs(fill = "Group", size = 12) +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour = "black")
    )  
}


############# Finding Clusters #############
cluster_finding <- function(raw_data, clu_ident){
  raw_data <- subset(pbmc.combined, idents = clu_ident)
  raw_data <- FindNeighbors(raw_data, dims = 1:21)
  raw_data <- FindClusters(raw_data, resolution = 0.5)
  raw_data <- RunUMAP(raw_data, dims = 1:21)
  
  #raw_data <- RunTSNE(raw_data, dims = 1:21)
  p1 <- DimPlot(raw_data, reduction = "umap", group.by = "orig.ident")
  p2 <- DimPlot(raw_data, reduction = "umap", label = T) 
  p3 <- DimPlot(raw_data, reduction = "umap", label = T, split.by = "orig.ident") + NoLegend()
  top_row <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
  p4 <- plot_grid(top_row, p3, labels = c('', 'C'), label_size = 12, ncol = 1)
  plot(p4)
  
  return(raw_data)
}

############################## Tumor Componants #####################
tumor_data <- cluster_finding(pbmc.combined, c(0, 1, 3, 5, 8, 9, 10, 11, 14, 20, 22, 7, 17))
tumor_meta <- tumor_data@meta.data %>% as.data.frame
component_stat(tumor_meta, 17)
tumor_markers <- FindAllMarkers(tumor_data, only.pos = 0.25, logfc.threshold = 0.25)
tumor_final_markers <- tumor_markers %>% group_by(cluster) %>% as.data.frame()

top10 <- tumor_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(tumor_data, features = top10$gene) + NoLegend()
################################ Microenvironment Clusters ###########################
micro_data <- cluster_finding(pbmc.combined, c(18, 15, 13, 16, 2, 12, 19, 6, 4, 21))
micro_meta <- micro_data@meta.data %>% as.data.frame
component_stat(micro_meta, 15)

micro_markers <- FindAllMarkers(micro_data, only.pos = 0.25, logfc.threshold = 0.25)
micro_final_markers <- micro_markers %>% group_by(cluster) %>% as.data.frame()
top10 <- micro_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(micro_data, features = top10$gene) + NoLegend()

dcast(micro_final_markers, cluster~gene) %>% View


Neu_Cells <- cluster_finding(pbmc.combined, 6)
Neu_meta <- B_Cells@meta.data %>% as.data.frame
component_stat(Neu_meta, 5)

Neu_markers <- FindAllMarkers(B_Cells, only.pos = 0.25, logfc.threshold = 0.25)
Neu_final_markers <- Neu_markers  %>% group_by(cluster) %>% as.data.frame()