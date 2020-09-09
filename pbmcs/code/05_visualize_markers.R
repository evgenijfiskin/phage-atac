library(Seurat)
library(BuenColors)
library(viridis)
library(dplyr)
load("../../../phage_atac_large_data_files/output/pbmcs/28July2020_analyzed_seurat_objects.rda")
#DimPlot(coembed, split.by = "orig.ident", label = TRUE)
#FindMarkers(coembed,ident.1 = "11")


make_feature_plot_chimera <- function(feature_name){
  
  DefaultAssay(coembed) <- "RNA"
  plot2 <- FeaturePlot(coembed, features = c(feature_name),
                       min.cutoff = "q01", max.cutoff = "q99",
                       split.by = "orig.ident", pt.size = 0.1)
  rna <- plot2[[2]] + theme_void() + theme(legend.position = "none")+ scale_color_viridis() + ggtitle(NULL)
  
  DefaultAssay(coembed) <- "ACTIVITY"
  plot2a <- FeaturePlot(coembed, features = c(feature_name),
                       min.cutoff = "q01", max.cutoff = "q90",
                       split.by = "orig.ident", pt.size = 0.1)
  atac <- plot2a[[1]] + theme_void() + theme(legend.position = "none")+ scale_color_viridis() + ggtitle(NULL)
  
  cowplot::ggsave2(cowplot::plot_grid(rna,atac, ncol = 1, scale = 0.9), file = paste0("../plots/markers/", feature_name, "_chimera2.pdf"), 
                   dpi = 500, width = 3, height = 6)
  feature_name
}

markers <- c("KLRD1", "CD8A", "LEF1", "MS4A1", "MAFB", "IL3RA")

sapply(markers, function(x){
  quantile(coembed@assays$RNA@data[x,coembed@meta.data$tech == "RNA"],c(0.01, 0.99))
})

sapply(markers, function(x){
  round(quantile(coembed@assays$ACTIVITY@data[x,coembed@meta.data$tech == "atac"],c(0.01, 0.9)),2)
})



lapply(markers, make_feature_plot_chimera)


make_feature_plot_atac <- function(feature_name){
  plot2 <- FeaturePlot(coembed, features = c(feature_name),
                       min.cutoff = "q00", max.cutoff = "q90",
                       split.by = "orig.ident")
  atac <- plot2[[1]] + theme_void() + theme(legend.position = "none") + scale_color_viridis() + ggtitle(NULL)
  rna <- plot2[[2]] + theme_void() + theme(legend.position = "none")+ scale_color_viridis() + ggtitle(NULL)
  cowplot::ggsave2(cowplot::plot_grid(atac,rna, ncol = 1, scale = 0.9), file = paste0("../plots/markers/", feature_name, "_ATAC.png"), 
                   dpi = 500, width = 3, height = 6)
  feature_name
}

batch1 <-   c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A", "CD19", "CCL4", "CCL5",
              "IL3RA", "CCR3", "HDC") 
ef <- c("LEF1", "CCR3","IKZF2", "NKG7", "KLRD1")
ef2 <- c("CST1", "MAFB", "GCSAML", "CST3")
genes_plot <- intersect(
  ef2,
  rownames(coembed@assays$RNA@scale.data)
)
lapply(genes_plot, make_feature_plot_atac)


