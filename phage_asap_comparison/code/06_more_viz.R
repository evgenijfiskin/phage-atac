library(Signac)
library(Seurat)
library(BuenColors)
library(Matrix)
library(pheatmap)
library(viridis)
library(dplyr)
load("../../../phage_atac_large_data_files/output/asap_phage_comparison/15Feb2021_analyzed_seurat_objects.rda")

# First, viz gene activity scores
coembed$channel <- substr(rownames(coembed@meta.data), 18, 18 )
pbmc.phage <- subset(coembed, channel %in% c(2,3) & celltype %in% c("CD14_Mono", "CD16_Mono", "Mem_CD4T", "Naive_CD4T", "Naive_CD8T", "Effector_CD8T" ))
DefaultAssay(pbmc.phage) <- "peaks"
pbmc.phage <- RunTFIDF(pbmc.phage)
pbmc.phage <- FindTopFeatures(pbmc.phage, min.cutoff = 'q50')
pbmc.phage <- RunSVD(pbmc.phage)

pbmc.phage <- RunUMAP(object = pbmc.phage, reduction = 'lsi', dims = 2:30)
pbmc.phage <- FindNeighbors(object = pbmc.phage, reduction = 'lsi', dims = 2:30)
pbmc.phage <- FindClusters(object = pbmc.phage, verbose = FALSE, algorithm = 3, resolution = 0.5)
DimPlot(pbmc.phage, label = TRUE)

FeaturePlot(pbmc.phage, c("NB17", "NB25"), max.cutoff = "q90") & scale_color_viridis()

df <- data.frame(CD4counts = pbmc.phage@assays$PDT@counts["NB17",], 
                 CD4CLR = pbmc.phage@assays$PDT@data["NB17",], 
                 CD8counts = pbmc.phage@assays$PDT@counts["NB25",], 
                 CD8CLR = pbmc.phage@assays$PDT@data["NB25",], 
                 pbmc.phage@meta.data)
df$celltype3 <- ifelse(df$seurat_clusters %in% c("0", "6"), "Monocytes",
                       ifelse(df$seurat_clusters %in% c("1", "2"), "CD4", 
                              ifelse(df$seurat_clusters %in% c("3", "4", "5"), "CD8", "other")))


dffilt <- df %>% dplyr::filter(
  !(celltype3 == "CD4" & CD4CLR < CD8CLR) &
    !(celltype3 == "CD8" & CD4CLR > CD8CLR) &
    celltype3 != "other" &
    !(CD4CLR >1.5 & CD8CLR > 2 & celltype3 == "CD8")
)

ggplot(dffilt, aes(x = log10(CD4counts), y = log10(CD8counts), color = celltype3)) +
  geom_point() + 
  scale_color_manual(values = c("firebrick", "dodgerblue2", "forestgreen"))


p1 <- ggplot(dffilt, aes(x = CD4CLR, fill = celltype3)) +
  geom_density(aes(color = celltype3), fill = NA,adjust = 2 ) +
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("firebrick", "dodgerblue2", "forestgreen")) +
  scale_fill_manual(values = c("firebrick", "dodgerblue2", "forestgreen")) +
  ggtitle("NB17") +  labs(x = "CLR CD4 counts") +
  theme(plot.title = element_text(size = 8)) + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/density_CD4clr.pdf", width = 2, height = 2)
table(dffilt$celltype)
