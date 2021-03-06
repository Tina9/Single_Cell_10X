library(Seurat)
library(dplyr)
library(patchwork)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
options(scipen = 20)
############ Identifing Cluster Changing ##########
setwd("~/Dropbox/BRCA1-PARPi-10X/Latest Version/cluster_finding/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest Version/cluster_finding/aggr1_2_3_combined_data.RData")
MetaData <- aggr1_2_3_combined_data@meta.data
######### Using identity number to plot #############
plot_stat <- MetaData %>% 
  group_by(orig.ident) %>%
  dplyr::count(seurat_clusters) %>%
  as.data.frame()

nb.cols <- 22
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)
pdf("clusterCount.pdf", height = 9, width = 6)
p1 <- ggplot(data = plot_stat, aes(x = orig.ident, y = n, fill = seurat_clusters)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = mycolors) +
        theme_classic() +
        xlab("") +
        ylab("Cell Numbers") +
        theme(text = element_text(size = 20),
              axis.text.x = element_text(angle = 45, hjust = 1)
              )
plot(p1)
dev.off()

#### Plot using cluster percentage
plot_percent_stat <- plot_stat %>%
  group_by(orig.ident) %>%
  mutate(per = (n/sum(n) * 100))

pdf("clusterPercent.pdf", height = 9, width = 9)
p2 <- ggplot(data = plot_percent_stat, aes(x = seurat_clusters, y = per, fill = orig.ident)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("darkblue", "red")) +
        theme_classic() +
        xlab("Cluster") +
        ylab("Cluster Percent") +
        labs(fill = "Group", size = 18) +
        theme(text = element_text(size = 18),
              axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
              axis.text.y = element_text(colour = "black")
        )
plot(p2)
dev.off()

write.table(plot_percent_stat, file = "stat.txt", row.names = F, sep = "\t", quote = F)
