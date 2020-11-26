library(Seurat)
library(dplyr)
library(patchwork)
library(cowplot)
library(ggplot2)
library(SingleR)
options(scipen = 20)

read_data <- function(data_dir, sample_name){
  pbmc.data <- Read10X(data_dir)
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = sample_name,
                             min.cells = 3, min.features = 200)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt")
  Violine_plot_name <- paste(sample_name, "_VIn.pdf", sep = "")
  pdf(Violine_plot_name, height = 6, width = 12)
  p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot(p1)
  dev.off()
  return(pbmc)
}

### Use density plot to filter cells
density_plot <- function(plot_data) {
  p1 <- ggplot(plot_data, aes(x=nCount_RNA, color=orig.ident)) +
    geom_density() +
    theme_bw() +
    theme(legend.position = "none")
  p2 <- ggplot(plot_data, aes(x = nFeature_RNA, color = orig.ident)) + 
    geom_density() +
    theme_bw() +
    theme(legend.position = "none")
  p3 <- ggplot(plot_data, aes(x = percent.mt, color = orig.ident)) +
    geom_density() +
    theme_bw()
  pdf("filter_density_plot.pdf", height = 4, width = 12)
  plot(p1 + p2 + p3)
  dev.off()
}

### fitering and normalizing data
prepare_data <- function(pbmc, min_feature, max_feature, mt_percent){
  pbmc <- subset(pbmc, subset = nFeature_RNA > min_feature & 
                   nFeature_RNA < max_feature &
                   percent.mt < mt_percent)
  pbmc <- NormalizeData(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
}

### Integrating data
data_integration <- function(merge_list){
  ################# Perform Intergration ####################
  immune.anchors <- FindIntegrationAnchors(object.list = merge_list)
  immune.combined <- IntegrateData(anchorset = immune.anchors)
  
  ################ Perform an integrated analysis ###########
  DefaultAssay(immune.combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  immune.combined <- ScaleData(immune.combined, verbose = F)
  immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = F)
  
  # t-SNE and Clustering
  immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
  immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
  immune.combined <- FindClusters(immune.combined, resolution = 0.5)
  # Visualization
  p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident")
  p2 <- DimPlot(immune.combined, reduction = "umap", label = T) 
  p3 <- DimPlot(immune.combined, reduction = "umap", label = T, split.by = "orig.ident") + NoLegend()
  top_row <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
  p4 <- plot_grid(top_row, p3, labels = c('', 'C'), label_size = 12, ncol = 1)
  pdf("cluster_seeking.pdf", height = 12, width = 16)
  plot(p4)
  dev.off()
  
  return(immune.combined)
}

### Seeking for clusters
seeking_clusters <- function(immune.combined){
  DefaultAssay(immune.combined) <- "integrated"
  all_markers <- FindAllMarkers(immune.combined, only.pos = 0.25, logfc.threshold = 0.25)
  final_markers <- all_markers %>% group_by(cluster) %>% as.data.frame()
  write.table(final_markers, "marker_genes.txt", sep = "\t", quote = F, row.names = F)
  
  top10 <- all_markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_logFC)
  pdf("combined_markers_heatmap.pdf", height = 20, width = 20)
  p2 <- DoHeatmap(immune.combined, features = top10$gene) + NoLegend()
  plot(p2)
  dev.off()
  
  # n = 0
  # for(i in seq(1, length(top10$gene), by = 9)){
  #   pdf(paste("cluster", n, ".pdf", sep = ""), height = 20, width = 16)
  #   p1 <- VlnPlot(immune.combined, features = top10$gene[i:(i+8)], split.by = "orig.ident")
  #   plot(p1)
  #   dev.off()
  #   n = n + 1
  # }
}

#### Using SingleR to annote Clusters
ClusterFinding <- function(ExpressData, metadata){
  mouse_refer <- ImmGenData()
  annote_res <- SingleR(test = ExpressData, ref = mouse_refer, labels = mouse_refer$label.main)
  metadata$Immlabels <- as.data.frame(annote_res)$labels
  return(metadata)
}

# #############################################################
# ############ Analysis for HP4568 ############################
# #############################################################
# setwd("~/Dropbox/BRCA1-PARPi-10X/Latest Version/HP4568/")
# control_data_dir <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4568_1/"
# control_name <- "untreated"
# exp_data_dir <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/aggr2_3/"
# exp_name <- "treated"
# control_data <- read_data(control_data_dir, control_name)
# exp_data <- read_data(exp_data_dir, exp_name)
# plot_data <- rbind(control_data@meta.data, exp_data@meta.data)
# density_plot(plot_data)
# data_list <- list(control_data, exp_data)
# names(data_list) <- c("untreated", "treated")
# data_list <- lapply(X = data_list, FUN = prepare_data, min_feature = 200,
#                     max_feature = 7000, mt_percent = 15)
# HP4568_combined_data <- data_integration(data_list)
# seeking_clusters(HP4568_combined_data)
# 
# ### Get cell cluster
# metadata <- HP4568_combined_data@meta.data %>% as.data.frame
# ### Get Expression Values
# ExpressData <- GetAssayData(HP4568_combined_data[["integrated"]], slot = "data")
# 
# ### Annotate Cell Clusters
# MetaData <- ClusterFinding(ExpressData, metadata)
# 
# ### Keep data for next usuage ##############
# save(HP4568_combined_data, file = "HP4568_combined_data.RData")
# write.table(MetaData, file = "HP4568_combined_data_combined_metadata.txt", sep = "\t", quote = F)
# write.table(as.data.frame(ExpressData), file = "HP4568 _features.csv", sep = ",", quote = F)
# 
# #############################################################
# ############ Analysis for HP4313 ############################
# #############################################################
# setwd("~/Dropbox/BRCA1-PARPi-10X/Latest Version/HP4313/")
# control_data_dir <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4313Pr/"
# control_name <- "untreated"
# exp_data_dir <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4313Re/"
# exp_name <- "treated"
# control_data <- read_data(control_data_dir, control_name)
# exp_data <- read_data(exp_data_dir, exp_name)
# plot_data <- rbind(control_data@meta.data, exp_data@meta.data)
# density_plot(plot_data)
# data_list <- list(control_data, exp_data)
# names(data_list) <- c("untreated", "treated")
# data_list <- lapply(X = data_list, FUN = prepare_data, min_feature = 200,
#                     max_feature = 7000, mt_percent = 15)
# HP4313_combined_data <- data_integration(data_list)
# seeking_clusters(HP4313_combined_data)
# 
# ### Get cell cluster
# metadata <- HP4313_combined_data@meta.data %>% as.data.frame
# ### Get Expression Values
# ExpressData <- GetAssayData(HP4313_combined_data[["integrated"]], slot = "data")
# 
# ### Annotate Cell Clusters
# MetaData <- ClusterFinding(ExpressData, metadata)
# 
# ### Keep data for next usuage ##############
# save(HP4313_combined_data, file = "HP4313_combined_data.RData")
# write.table(MetaData, file = "HP4313_combined_metadata.txt", sep = "\t", quote = F)
# write.table(as.data.frame(ExpressData), file = "HP4313_features.csv", sep = ",", quote = F)
# 
# #############################################################
# ############ Analysis for HP4568 ############################
# #############################################################
# setwd("~/Dropbox/BRCA1-PARPi-10X/Latest Version/HP4568/")
# control_data_dir <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4568Pr/"
# control_name <- "untreated"
# exp_data_dir <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4568Re/"
# exp_name <- "treated"
# control_data <- read_data(control_data_dir, control_name)
# exp_data <- read_data(exp_data_dir, exp_name)
# plot_data <- rbind(control_data@meta.data, exp_data@meta.data)
# density_plot(plot_data)
# data_list <- list(control_data, exp_data)
# names(data_list) <- c("untreated", "treated")
# data_list <- lapply(X = data_list, FUN = prepare_data, min_feature = 200,
#                     max_feature = 7000, mt_percent = 15)
# HP4568_combined_data <- data_integration(data_list)
# seeking_clusters(HP4568_combined_data)
# 
# ### Get cell cluster
# metadata <- HP4568_combined_data@meta.data %>% as.data.frame
# ### Get Expression Values
# ExpressData <- GetAssayData(HP4568_combined_data[["integrated"]], slot = "data")
# 
# ### Annotate Cell Clusters
# MetaData <- ClusterFinding(ExpressData, metadata)
# 
# ### Keep data for next usuage ##############
# save(HP4568_combined_data, file = "HP4568_combined_data.RData")
# write.table(MetaData, file = "HP4568_combined_metadata.txt", sep = "\t", quote = F)
# write.table(as.data.frame(ExpressData), file = "HP4568_features.csv", sep = ",", quote = F)


HP4540_Pre_data_dir <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4540_1/"
HP4540_Pre_name <- "HP4540Pre"
HP4540_Pre_data <- read_data(HP4540_Pre_data_dir, HP4540_Pre_name)
HP4540_Tr_data_dir <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/aggr2_3/"
HP4540_Tr_name <- "HP4540Tr"
HP4540_Tr_data <- read_data(HP4540_Tr_data_dir, HP4540_Tr_name)
HP4313_Pre_data_dir <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4313Pr/"
HP4313_Pre_name <- "HP4313Pre"
HP4313_Pre_data <- read_data(HP4313_Pre_data_dir, HP4313_Pre_name)
HP4313_Tr_data_dir <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4313Re/"
HP4313_Tr_name <- "HP4313Tr"
HP4313_Tr_data <- read_data(HP4313_Tr_data_dir, HP4313_Tr_name)
HP4568_Pre_data_dir <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4568Pr/"
HP4568_Pre_name <- "HP4568Pre"
HP4568_Pre_data <- read_data(HP4568_Pre_data_dir, HP4568_Pre_name)
HP4568_Tr_data_dir <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4568Re/"
HP4568_Tr_name <- "HP4568Tr"
HP4568_Tr_data <- read_data(HP4568_Tr_data_dir, HP4568_Tr_name)
plot_data <- rbind(HP4540_Pre_data@meta.data,
                   HP4540_Tr_data@meta.data,
                   HP4313_Pre_data@meta.data,
                   HP4313_Tr_data@meta.data,
                   HP4568_Pre_data@meta.data,
                   HP4568_Tr_data@meta.data)
density_plot(plot_data)
data_list <- list(HP4540_Pre_data,
                  HP4540_Tr_data,
                  HP4313_Pre_data,
                  HP4313_Tr_data,
                  HP4568_Pre_data,
                  HP4568_Tr_data)
names(data_list) <- c("HP4540_Pre",
                      "HP4540_Tr",
                      "HP4313_Pre",
                      "HP4313_Tr",
                      "HP4568_Pre",
                      "HP4568_Tr")
data_list <- lapply(X = data_list, FUN = prepare_data, min_feature = 200,
                    max_feature = 7000, mt_percent = 15)
total_combined_data <- data_integration(data_list)
seeking_clusters(total_combined_data)

### Get cell cluster
metadata <- total_combined_data@meta.data %>% as.data.frame
### Get Expression Values
ExpressData <- GetAssayData(total_combined_data[["integrated"]], slot = "data")

### Annotate Cell Clusters
MetaData <- ClusterFinding(ExpressData, metadata)

### Keep data for next usuage ##############
save(total_combined_data, file = "total_combined_data.RData")
write.table(MetaData, file = "total_combined_data_combined_metadata.txt", sep = "\t", quote = F)
write.table(as.data.frame(ExpressData), file = "total _features.csv", sep = ",", quote = F)

