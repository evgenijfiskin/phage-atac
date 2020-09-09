library(Seurat)
library(BuenColors)
library(dplyr)

load("../../../phage_atac_large_data_files/output/pbmcs/28July2020_analyzed_seurat_objects.rda")

coembed@meta.data$pct_in_peaks <- coembed@meta.data$peak_region_fragments/coembed@meta.data$passed_filters
DimPlot(coembed, group.by = "seurat_clusters", split.by = "tech", label = TRUE)

FeaturePlot(pbmc.cite, features = c("adt_CD16", "adt_CD4", "adt_CD8"),
            min.cutoff = "q01", max.cutoff = "q99", ncol = 3, reduction = "umapd")

Idents(pbmc.atac) <-  coembed@meta.data$seurat_clusters[coembed@meta.data$tech == "atac"]
Idents(pbmc.cite) <-  coembed@meta.data$seurat_clusters[coembed@meta.data$tech == "RNA"]
dim(pbmc.cite)
dim(pbmc.atac)

#FindMarkers(pbmc.atac, ident.1 = "11", assay = "ACTIVITY") %>% head(50)
FindMarkers(pbmc.cite, ident.1 = "12") %>% head(10)
#FeaturePlot(coembed, features = c("S100A4"))

# rename
vec <- c(
  "0" = "CD14_Mono",
  "1" = "Naive_CD4T",
  "2" = "Memory_CD4T",
  "3" = "NK_cell",
  "4" = "Bcell",
  "5" = "Effector_CD8T",
  "6" = "CD16_Mono",
  "7" = "Naive_CD8T",
  "8" = "Cytotoxic_CD8T",
  "9" = "mDC",
  "10" = "Plasma_cells",
  "11" = "Granulocyte",
  "12" = "pDC",
  "13" = "Megakaryocyte"
)
coembed$celltype <- vec[as.character(coembed$seurat_clusters)]
pbmc.atac$celltype <- vec[as.character(Idents(pbmc.atac))]
pbmc.cite$celltype <- vec[as.character(Idents(pbmc.cite))]

pU1 <- DimPlot(pbmc.cite, reduction = "umapd", group.by = "celltype", pt.size = 0.3) +
  theme_void() + 
  scale_color_manual(values = c(jdb_palette("corona")[1:14])) +
  theme(legend.position = "none") 

pU2 <- DimPlot(pbmc.atac, reduction = "umapd", group.by = "celltype", pt.size = 0.3) +
  theme_void() +
  scale_color_manual(values = c(jdb_palette("corona")[1:14])) +
  theme(legend.position = "none") 

cowplot::ggsave2(cowplot::plot_grid(pU1, pU2, ncol =1), 
                 file = "../plots/UMAP_clusters_stack.png", width = 3, height = 6, dpi = 500)

atac_df <- data.frame(
  cluster = vec[as.character(Idents(pbmc.atac))],
  t(pbmc.atac@assays$ADT@scale.data)
) %>% group_by(cluster) %>% summarize_all(mean)

rna_df <- data.frame(
  cluster = vec[as.character(Idents(pbmc.cite))],
  t(pbmc.cite@assays$ADT@scale.data)
) %>% group_by(cluster) %>% summarize_all(mean)

mdf <- merge(atac_df, rna_df, by = "cluster")
color_vec <- jdb_palette("corona", 14); names(color_vec) <- mdf[[1]]
ps0 <- ggplot(mdf, aes(x = CD16.x, y = CD16.y, color = cluster)) +
  geom_point(size = 0.9) + pretty_plot(fontsize = 8) + L_border() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = color_vec) 
cowplot::ggsave2(g_legend(ps0), file = "../plots/pbmc_legend.pdf", 
                 width = 2, height = 5)

mdf <- mdf %>% dplyr::filter(!(cluster %in% c("Megakaryocyte", "Granulocyte")))

min(mdf[,c(4,7)])
max(mdf[,c(4,7)])
cor.test(mdf[,c(4)],mdf[,c(7)])

ps1 <- ggplot(mdf, aes(x = CD16.x, y = CD16.y, color = cluster)) +
  geom_point(size = 0.9) + pretty_plot(fontsize = 8) + L_border() +
  labs(x = "CD16 Phage", y = "CD16 CITE") +
  scale_color_manual(values = color_vec) +
  scale_x_continuous(breaks = c(0,1), limits = c(-0.62, 1.75)) +
  scale_y_continuous(breaks = c(0,1), limits = c(-0.62, 1.75)) +
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

min(mdf[,c(2,5)])
max(mdf[,c(2,5)])
cor.test(mdf[,c(2)],mdf[,c(5)])


ps2 <- ggplot(mdf, aes(x = CD4.x, y = CD4.y, color = cluster)) +
  geom_point(size = 0.9) + pretty_plot(fontsize = 8) + L_border() +
  labs(x = "CD4 Phage", y = "CD4 CITE") +
  scale_color_manual(values = color_vec) +
  scale_x_continuous(breaks = c(0,1), limits = c(-1.0, 1.05)) +
  scale_y_continuous(breaks = c(0,1), limits = c(-1.0, 1.05)) + theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

min(mdf[,c(3,6)])
max(mdf[,c(3,6)])
cor.test(mdf[,c(3)],mdf[,c(6)])

ps3 <- ggplot(mdf, aes(x = CD8.x, y = CD8.y, color = cluster)) +
  geom_point(size = 0.9) + pretty_plot(fontsize = 8) + L_border() +
  labs(x = "CD8 Phage", y = "CD8 CITE") +
  scale_color_manual(values = color_vec) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0,1,2), limits = c(-0.66, 2.15)) +
  scale_y_continuous(breaks = c(0,1,2), limits = c(-0.66, 2.15)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

cowplot::ggsave2(cowplot::plot_grid(ps1, ps2, ps3, ncol =1), file = "../plots/scatter_ADTs.pdf", width = 1, height = 3)

# Now look in UMAP space

pP1 <- FeaturePlot(pbmc.atac, features = c("adt_CD16"),
                   min.cutoff = "q01", max.cutoff = "q99", reduction = "umapd", pt.size = 0.3) + theme_void() +
  ggtitle("") + theme(plot.title = element_text(size = 8)) +
  theme(legend.position = "none") + scale_color_viridis()

pP2 <- FeaturePlot(pbmc.atac, features = c("adt_CD4"),
                   min.cutoff = "q01", max.cutoff = "q99", reduction = "umapd", pt.size = 0.3) + theme_void() +
  ggtitle("") + theme(plot.title = element_text(size = 8)) +
  theme(legend.position = "none") + scale_color_viridis()

pP3 <- FeaturePlot(pbmc.atac, features = c("adt_CD8"),
                   min.cutoff = "q01", max.cutoff = "q99", reduction = "umapd", pt.size = 0.3) + theme_void() +
  ggtitle("") + theme(plot.title = element_text(size = 8)) +
  theme(legend.position = "none") + scale_color_viridis()

pC1 <- FeaturePlot(pbmc.cite, features = c("adt_CD16"),
                   min.cutoff = "q01", max.cutoff = "q99", reduction = "umapd", pt.size = 0.3) + theme_void() +
  ggtitle("") + theme(plot.title = element_text(size = 8)) +
  theme(legend.position = "none")  + scale_color_viridis()

pC2 <- FeaturePlot(pbmc.cite, features = c("adt_CD4"),
                   min.cutoff = "q01", max.cutoff = "q99", reduction = "umapd", pt.size = 0.3) + theme_void() +
  ggtitle("") + theme(plot.title = element_text(size = 8)) +
  theme(legend.position = "none") + scale_color_viridis()

pC3 <- FeaturePlot(pbmc.cite, features = c("adt_CD8"),
                   min.cutoff = "q01", max.cutoff = "q99", reduction = "umapd", pt.size = 0.3) + theme_void() +
  ggtitle("") + theme(plot.title = element_text(size = 8)) +
  theme(legend.position = "none") + scale_color_viridis()

cowplot::ggsave2(cowplot::plot_grid(pC1, pC2, pC3, pP1, pP2, pP3, ncol =3), 
                 file = "../plots/UMAP_ADTs_stack.png", width = 9, height = 6, dpi = 500)


sapply(c("CD8", "CD4", "CD16"), function(x){
  quantile(pbmc.atac@assays$ADT@data[x,],c(0.01, 0.99))
})
