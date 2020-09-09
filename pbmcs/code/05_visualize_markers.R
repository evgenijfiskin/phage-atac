library(Seurat)
library(BuenColors)
library(viridis)
library(dplyr)
load("../../../phage_atac_large_data_files/output/pbmcs/28July2020_analyzed_seurat_objects.rda")
#DimPlot(coembed, split.by = "orig.ident", label = TRUE)
#FindMarkers(coembed,ident.1 = "11")

DefaultAssay(coembed) <- "RNA"

make_feature_plot <- function(feature_name){
  plot2 <- FeaturePlot(coembed, features = c(feature_name),
                       min.cutoff = "q01", max.cutoff = "q99",
                       split.by = "orig.ident")
  atac <- plot2[[1]] + theme_void() + theme(legend.position = "none") + scale_color_viridis() + ggtitle(NULL)
  rna <- plot2[[2]] + theme_void() + theme(legend.position = "none")+ scale_color_viridis() + ggtitle(NULL)
  cowplot::ggsave2(cowplot::plot_grid(atac,rna, ncol = 1, scale = 0.9), file = paste0("../plots/markers/", feature_name, ".png"), 
                   dpi = 500, width = 3, height = 6)
  feature_name
}

batch1 <-   c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A", "CD19", "CCL4", "CCL5",
              "IL3RA", "CCR3", "HDC") 
ef <- c("LEF1", "CCR3","IKZF2", "NKG7", "KLRD1")
genes_plot <- intersect(
  ef,
  rownames(coembed@assays$RNA@scale.data)
)
lapply(genes_plot, make_feature_plot)

