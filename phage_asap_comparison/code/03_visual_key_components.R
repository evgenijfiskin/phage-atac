library(Signac)
library(Seurat)
library(BuenColors)
library(Matrix)
library(pheatmap)
library(viridis)
library(dplyr)
load("../../../phage_atac_large_data_files/output/asap_phage_comparison/15Feb2021_analyzed_seurat_objects.rda")


# First, viz gene activity scores
pbmc.phage.asap <- subset(pbmc.atac, channel == 2)
make_feature_plot_chimera <- function(feature_name){
  
  DefaultAssay(coembed) <- "RNA"
  plot2 <- FeaturePlot(coembed, features = c(feature_name),
                       min.cutoff = "q01", max.cutoff = "q99",
                       split.by = "orig.ident", pt.size = 0.1)
  rna <- plot2[[2]] + theme_void() + theme(legend.position = "none") + ggtitle(NULL)
  
  DefaultAssay(coembed) <- "ACTIVITY"
  plot2a <- FeaturePlot(pbmc.phage.asap, features = c(feature_name),
                        min.cutoff = "q01", max.cutoff = "q90", reduction = "umapd",
                        split.by = "orig.ident", pt.size = 0.1)
  atac <- plot2a[[1]] + theme_void() + theme(legend.position = "none")+  ggtitle(NULL)
  
  cowplot::ggsave2(cowplot::plot_grid(rna,atac, ncol = 1, scale = 0.9), file = paste0("../plots/markers/", feature_name, "_chimera2-update24feb.pdf"), 
                  width = 3, height = 6)
  feature_name
}
make_feature_plot_chimera("CD8A")
markers <- c("KLRD1", "CD8A", "LEF1", "MS4A1", "MAFB", "IL3RA", "CD4", "CEACAM4", "FCGR3A")
lapply(markers, make_feature_plot_chimera)

sapply(markers, function(x){
  quantile(coembed@assays$RNA@data[x,coembed@meta.data$tech == "RNA"],c(0.01, 0.99))
})

sapply(markers, function(x){
  round(quantile(coembed@assays$ACTIVITY@data[x,coembed@meta.data$tech4 == "ASAPandPHAGE"],c(0.01, 0.9)),2)
})

make_feature_plot_chimera("FCGR3A")



ct <- sort(unique(coembed$celltype))[1:10] # rm Platelets and Progenitors
cite_mat <- sapply(unique(ct), function(celltype){
  celltype_vec <- coembed$celltype == celltype
  rowMeans(coembed@assays$ADT@data[,celltype_vec & coembed$tech4 == "CITE"])
})

phage_mat <- sapply(unique(ct), function(celltype){
  celltype_vec <- coembed$celltype == celltype
  rowMeans(coembed@assays$PDT@data[c("NB17", "NB18", "NB25"),celltype_vec & coembed$tech4 == "PHAGE"])
})

asap_mat <- sapply(unique(ct), function(celltype){
  celltype_vec <- coembed$celltype == celltype
  rowMeans(coembed@assays$ASAP@data[,celltype_vec & coembed$tech4 == "ASAP"])
})

cophage_mat <- sapply(unique(ct), function(celltype){
  celltype_vec <- coembed$celltype == celltype
  rowMeans(coembed@assays$PDT@data[c("NB17", "NB18", "NB25"),celltype_vec & coembed$tech4 == "ASAPandPHAGE"])
})
coasap_mat <- sapply(unique(ct), function(celltype){
  celltype_vec <- coembed$celltype == celltype
  rowMeans(coembed@assays$ASAP@data[,celltype_vec & coembed$tech4 == "ASAPandPHAGE"])
})
rownames(cophage_mat) <- paste0("co", rownames(cophage_mat))
mat <- cor(cbind(t(phage_mat), t(asap_mat), t(cite_mat), t(cophage_mat)), use =  "pairwise.complete")

matco <- cor(cbind(t(coembed@assays$PDT@data[c("NB17", "NB18", "NB25", "NB112"),coembed$tech4 == "ASAPandPHAGE"]),
                   t(coembed@assays$ASAP@data[,coembed$tech4 == "ASAPandPHAGE"])))

mat_one_channel <- c("NB17", "NB18", "NB25", "NB112", "asapCD4", "asapCD19", "asapCD11c", "asapCD3", "asapCD14")
pp <- pheatmap(matco[mat_one_channel,mat_one_channel], cluster_rows = TRUE, cluster_cols = TRUE,
               color = viridis(500), fontsize = 6)
pp
pdf("../plots/protein_modalities_cluster_correlation_big-ceacam.pdf", width = 5.5, height = 5)
pp
dev.off()
order <- c("NB25", "asapCD8", "coNB25", "CD8", "NB18", "asapCD16", "coNB18", "CD16",
           "NB17", "asapCD4", "coNB17","CD4")
pp <- pheatmap(mat[order,order], cluster_rows = TRUE, cluster_cols = TRUE,
               color = viridis(500), fontsize = 6)
pp
pdf("../plots/protein_modalities_cluster_correlation_cosmall-clsuter.pdf", width = 3.5, height = 3)
pp
dev.off()

Idents(coembed) <- coembed$celltype

lapply(c("asapCD8", "CD8", "NB17", "NB18", "NB25", "asapCD4", "CD4", "CD16", "asapCD16"), function(feat){
  p1 <- VlnPlot(object = coembed, features = c(feat), split.by = 'tech4', split.plot = TRUE, pt.size = 0,)
  cowplot::ggsave2(p1, file = paste0("../plots/EF_15Feb/violin_", feat, "_15Feb2021.pdf"), width = 6, height = 4)
})


p1 <- DimPlot(coembed, group.by = "celltype", pt.size = 0.3, split.by = "tech4") +
  theme_void() + 
  scale_color_manual(values = c(jdb_palette("corona")[1:14])) +
  theme(legend.position = "none") 

cowplot::ggsave2(p1, file = "../plots/UMAP_celltypes_row.png", width = 11, height = 3, dpi = 500)
cowplot::ggsave2(p1, file = "../plots/UMAP_celltypes_row.pdf", width = 11, height = 3)

legend <- g_legend(DimPlot(coembed, group.by = "celltype", pt.size = 0.3, split.by = "tech4") + scale_color_manual(values = c(jdb_palette("corona")[1:14])))

cowplot::ggsave2(legend, file = "../plots/color_legend.pdf", height = 4, width = 2)

library(viridis)

process_protein_plot <- function(feature, tech_filt, assay){
  viz_obj <- subset(coembed, tech4 == tech_filt)
  DefaultAssay(viz_obj) <- assay
  viz_obj <- suppressMessages(suppressWarnings(NormalizeData(viz_obj, assay = assay, normalization.method = "CLR", verbose = FALSE)))
  p1 <- FeaturePlot(viz_obj, features = c(feature), pt.size = 0.3, max.cutoff = "q99", min.cutoff = "q01") + theme_void() +
    ggtitle("") + theme(plot.title = element_text(size = 8))  + scale_color_viridis()
  cowplot::ggsave2(p1, file = paste0("../plots/EF_15Feb/UMAP_", feature, "_", tech_filt, ".pdf"), width = 3.6, height = 3)
  print(round(quantile(viz_obj@assays[[viz_obj@active.assay]]@data[feature,], probs = c(0.01, 0.99)),2))
}

process_protein_plot("CD4", "CITE", "ADT")
process_protein_plot("CD8", "CITE", "ADT")
process_protein_plot("CD16", "CITE", "ADT")

process_protein_plot("asapCD4", "ASAP", "ASAP")
process_protein_plot("asapCD4", "ASAPandPHAGE", "ASAP")
process_protein_plot("asapCD8", "ASAP", "ASAP")
process_protein_plot("asapCD8", "ASAPandPHAGE", "ASAP")
process_protein_plot("asapCD16", "ASAP", "ASAP")
process_protein_plot("asapCD16", "ASAPandPHAGE", "ASAP")

process_protein_plot("NB17", "PHAGE", "PDT")
process_protein_plot("NB17", "ASAPandPHAGE", "PDT")
process_protein_plot("NB18", "PHAGE", "PDT")
process_protein_plot("NB18", "ASAPandPHAGE", "PDT")
process_protein_plot("NB25", "PHAGE", "PDT")
process_protein_plot("NB25", "ASAPandPHAGE", "PDT")

process_protein_plot("NB112", "PHAGE", "PDT")
process_protein_plot("NB112", "ASAPandPHAGE", "PDT")


# Others
process_protein_plot("asapCD45", "ASAP", "ASAP")
process_protein_plot("asapCD45", "ASAPandPHAGE", "ASAP")

process_protein_plot("asapCD56", "ASAP", "ASAP")
process_protein_plot("asapCD56", "ASAPandPHAGE", "ASAP")

process_protein_plot("asapCD19", "ASAP", "ASAP")
process_protein_plot("asapCD19", "ASAPandPHAGE", "ASAP")

process_protein_plot("asapCD11c", "ASAP", "ASAP")
process_protein_plot("asapCD11c", "ASAPandPHAGE", "ASAP")

process_protein_plot("asapCD14", "ASAP", "ASAP")
process_protein_plot("asapCD14", "ASAPandPHAGE", "ASAP")

process_protein_plot("asapCD3", "ASAP", "ASAP")
process_protein_plot("asapCD3", "ASAPandPHAGE", "ASAP")


# Make ecdf plots

pa_df <- data.frame(
  CD4 = coembed@assays$PDT@counts["NB17", substr(rownames(coembed@meta.data), 18, 18) == "3"],
  CD8 = coembed@assays$PDT@counts["NB17", substr(rownames(coembed@meta.data), 18, 18) == "3"],
  CD16 = coembed@assays$PDT@counts["NB18", substr(rownames(coembed@meta.data), 18, 18) == "3"],
  cl =  coembed@meta.data$celltype[ substr(rownames(coembed@meta.data), 18, 18) == "3"]
)

# Make ecdf plots
pecdf <- ggplot(pa_df %>% filter(cl %in% c("CD14_Mono", "CD16_Mono")), aes(x = CD4 + 1)) + 
  stat_ecdf(color = "firebrick") + scale_x_log10(limits = c(1,500)) +
  labs(x = "log10 CD4 + 1 counts", y = "Cumulative Distribution")+ 
  pretty_plot(fontsize = 8) + L_border()
cowplot::ggsave2(pecdf, filename = "../plots/ecdf_Mono_CD4_counts_newphage_filter.pdf", width = 2.2, height = 2.2)

pecdf <- ggplot(pa_df %>% filter(cl %in% c("Mem_CD4T", "Naive_CD4T")), aes(x = CD4 + 1)) + 
  stat_ecdf(color = "firebrick") + scale_x_log10(limits = c(1,500)) +
  labs(x = "log10 CD4 + 1 counts", y = "Cumulative Distribution")+
  pretty_plot(fontsize = 8) + L_border()
cowplot::ggsave2(pecdf, filename = "../plots/ecdf_CD4T_CD4_counts_newphage_filter.pdf", width = 2.2, height = 2.2)


pecdf <- ggplot(pa_df %>% filter(cl %in% c("NK_cell")), aes(x = CD16 + 1)) + 
  stat_ecdf(color = "firebrick") + scale_x_log10(limits = c(1,500)) +
  labs(x = "log10 CD16 + 1 counts", y = "Cumulative Distribution")+ 
  pretty_plot(fontsize = 8) + L_border()
cowplot::ggsave2(pecdf, filename = "../plots/ecdf_NK_CD16_counts_newphage_filter.pdf", width = 2.2, height = 2.2)

pecdf <- ggplot(pa_df %>% filter(cl %in% c("CD16_Mono")), aes(x = CD16 + 1)) + 
  stat_ecdf(color = "firebrick") + scale_x_log10(limits = c(1,500)) +
  labs(x = "log10 CD16 + 1 counts", y = "Cumulative Distribution")+
  pretty_plot(fontsize = 8) + L_border()
cowplot::ggsave2(pecdf, filename = "../plots/ecdf_CD16Mono_CD16_counts_newphage_filter.pdf", width = 2.2, height = 2.2)



ct <- sort(unique(coembed$celltype))[1:10] # rm Platelets and Progenitors
cophage_mat <- sapply(unique(ct), function(celltype){
  cells <- colnames(coembed)[(coembed$celltype == celltype & coembed$tech4 == "ASAPandPHAGE")]
  rowMeans(coembed@assays$PDT@data[c("NB17", "NB18", "NB25", "NB112"),cells ])
})
coasap_mat <- sapply(unique(ct), function(celltype){
  celltype_vec <- coembed$celltype == celltype
  rowMeans(coembed@assays$ASAP@data[,celltype_vec & coembed$tech4 == "ASAPandPHAGE"])
})

keep <- c("asapCD3", "asapCD4", "asapCD14", "asapCD11c", "asapCD19", "NB17", "NB25", "NB18", "NB112")
pp <- pheatmap(cor(cbind(t(cophage_mat), t(coasap_mat))[,keep]), cluster_rows = TRUE, cluster_cols = TRUE,
               color = viridis(500), fontsize = 6)
pp
pdf("../plots/all_phage_asap_clusters-keep.pdf", width = 5.5, height = 5)
pp
dev.off()
