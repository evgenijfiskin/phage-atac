library(Seurat)
library(BuenColors)
library(dplyr)

load("../../../phage_atac_large_data_files/output/pbmcs/28July2020_analyzed_seurat_objects.rda")

FeaturePlot(pbmc.cite, features = c("adt_CD16", "adt_CD4", "adt_CD8"),
            min.cutoff = "q01", max.cutoff = "q99", ncol = 3, reduction = "umapd")


FindMarkers(pbmc.cite, ident.1 = "4", ident.2 = "7") %>% head(50)
FeaturePlot(coembed, features = c("CTLA4"))

# rename
vec <- c(
  "0" = "Naive_CD4",
  "1" = "CD14_Mono",
  "2" = "NK_cells",
  "3" = "Memory_CD4",
  "4" = "Effector_CD8",
  "5" = "Naive_Bcell",
  "6" = "Memory_Bcell",
  "7" = "Cytotoxic_Tcell",
  "8" = "Naive_CD8",
  "9" = "CD16_mono",
  "10" = "CD56bright_NK",
  "11" = "mDC",
  "12" = "Treg",
  "13" = "Platelet"
  )
coembed$celltype <- vec[as.character(coembed$seurat_clusters_m)]
pbmc.atac$celltype <- vec[as.character(pbmc.atac$predicted.id)]
pbmc.cite$celltype <- vec[as.character(pbmc.cite$seurat_clusters)]

pU1 <- DimPlot(pbmc.cite, reduction = "umapd", group.by = "celltype") +
  theme_void() + theme(plot.title = element_text(size = 8)) +
  scale_color_manual(values = c(jdb_palette("corona")[1:13], "lightgrey")) +
  theme(legend.position = "none") + ggtitle("CITE-seq (Stoeckius 2017)")

pU2 <- DimPlot(pbmc.atac, reduction = "umapd", group.by = "celltype") +
  theme_void() + theme(plot.title = element_text(size = 8)) +
  scale_color_manual(values = c(jdb_palette("corona")[1:13], "lightgrey")) +
  theme(legend.position = "none") + ggtitle("Phage ATAC") 

cowplot::ggsave2(cowplot::plot_grid(pU1, pU2, ncol =1), 
                 file = "../plots/UMAP_clusters_stack.pdf", width = 1.5, height = 3)

atac_df <- data.frame(
  cluster = vec[as.character(pbmc.atac$predicted.id)],
  t(pbmc.atac@assays$ADT@scale.data)
) %>% group_by(cluster) %>% summarize_all(mean)

rna_df <- data.frame(
  cluster = vec[as.character(pbmc.cite$seurat_clusters)],
  t(pbmc.cite@assays$ADT@scale.data)
) %>% group_by(cluster) %>% summarize_all(mean)

mdf <- merge(atac_df, rna_df, by = "cluster")

ps1 <- ggplot(mdf, aes(x = CD16.x, y = CD16.y, color = cluster)) +
  geom_point(size = 0.9) + pretty_plot(fontsize = 8) + L_border() +
  labs(x = "CD16 Phage ATAC", y = "CD16 CITE") +
  scale_color_manual(values = jdb_palette("corona")) +
  scale_x_continuous(breaks = c(0,1,2), limits = c(-0.5, 2)) +
  scale_y_continuous(breaks = c(0,1,2), limits = c(-0.5, 2)) +
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

ps2 <- ggplot(mdf, aes(x = CD4.x, y = CD4.y, color = cluster)) +
  geom_point(size = 0.9) + pretty_plot(fontsize = 8) + L_border() +
  labs(x = "CD4 Phage ATAC", y = "CD4 CITE") +
  scale_color_manual(values = jdb_palette("corona")) +
  scale_x_continuous(breaks = c(0,1,2), limits = c(-1.05, 2)) +
  scale_y_continuous(breaks = c(0,1,2), limits = c(-1.05, 2)) + theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

ps3 <- ggplot(mdf, aes(x = CD8.x, y = CD8.y, color = cluster)) +
  geom_point(size = 0.9) + pretty_plot(fontsize = 8) + L_border() +
  labs(x = "CD8 Phage ATAC", y = "CD8 CITE") +
  scale_color_manual(values = jdb_palette("corona")) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0,1,2), limits = c(-0.65, 2.5)) +
  scale_y_continuous(breaks = c(0,1,2), limits = c(-0.65, 2.5)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

cowplot::ggsave2(cowplot::plot_grid(ps1, ps2, ps3, ncol =1), file = "../plots/scatter_ADTs.pdf", width = 1.1, height = 3)

pP1 <- FeaturePlot(pbmc.atac, features = c("adt_CD16"),
            min.cutoff = "q01", max.cutoff = "q99", reduction = "umapd") + theme_void() +
  ggtitle("Phage CD16") + theme(plot.title = element_text(size = 8)) +
  theme(legend.position = "none")

pP2 <- FeaturePlot(pbmc.atac, features = c("adt_CD4"),
                   min.cutoff = "q01", max.cutoff = "q99", reduction = "umapd") + theme_void() +
  ggtitle("Phage CD4") + theme(plot.title = element_text(size = 8)) +
  theme(legend.position = "none")

pP3 <- FeaturePlot(pbmc.atac, features = c("adt_CD8"),
                   min.cutoff = "q01", max.cutoff = "q99", reduction = "umapd") + theme_void() +
  ggtitle("Phage CD8") + theme(plot.title = element_text(size = 8)) +
  theme(legend.position = "none")

pC1 <- FeaturePlot(pbmc.cite, features = c("adt_CD16"),
                   min.cutoff = "q01", max.cutoff = "q99", reduction = "umapd") + theme_void() +
  ggtitle("CITE CD16") + theme(plot.title = element_text(size = 8)) +
  theme(legend.position = "none")

pC2 <- FeaturePlot(pbmc.cite, features = c("adt_CD4"),
                   min.cutoff = "q01", max.cutoff = "q99", reduction = "umapd") + theme_void() +
  ggtitle("CITE CD4") + theme(plot.title = element_text(size = 8)) +
  theme(legend.position = "none")

pC3 <- FeaturePlot(pbmc.cite, features = c("adt_CD8"),
                   min.cutoff = "q01", max.cutoff = "q99", reduction = "umapd") + theme_void() +
  ggtitle("CITE CD8") + theme(plot.title = element_text(size = 8)) +
  theme(legend.position = "none")

cowplot::ggsave2(cowplot::plot_grid(pC1, pC2, pC3, pP1, pP2, pP3, ncol =3), 
                 file = "../plots/UMAP_ADTs_stack.pdf", width = 4.5, height = 3)

cowplot::ggsave2(g_legend(ps1), file = "../plots/pbmc_legend.pdf", 
                 width = 2, height = 5)
