library(Signac)
library(Seurat)
library(BuenColors)
library(Matrix)
library(pheatmap)
library(viridis)
library(dplyr)
# Import data and subset to interesting variants from previous scripts
load("../../../phage_atac_large_data_files/output/asap_phage_comparison/15Feb2021_analyzed_seurat_objects.rda")

col_vec <- jdb_palette("corona")[1:12]
names(col_vec) <- sort(unique(unname(c(
  "0" = "Naive_CD4T",
  "1" = "Mem_CD4T",
  "2" = "CD14_Mono",
  "3" = "Naive_Bcell",
  "4" = "NK_cell",
  "5" = "Effector_CD8T",
  "6" = "Naive_CD8T",
  "7" = "CD14_Mono",
  "8" = "Naive_CD4T",
  "9" = "Memory_Bcell",
  "10" = "CD16_Mono",
  "11" = "DCs",
  "12" = "Progenitors",
  "13" = "Platlets"
))))

# Make a dedicated version for CITE-seq that removes the weird doublet cluster

plot_me <- subset(coembed, tech4 =="CITE" & seurat_clusters != "7" & seurat_clusters != "13" & seurat_clusters != "12" )
DefaultAssay(plot_me) <- assay
Idents(plot_me) <- plot_me$celltype
plot_me <- NormalizeData(plot_me, assay = assay, normalization.method = "CLR")

pgo <- RidgePlot(plot_me, features  = "CD4", sort = "increasing") + 
  ggtitle(NULL) +   scale_fill_manual(values = col_vec) +
  theme(legend.position = "none") 

cowplot::ggsave2(pgo,
                 file = paste0("../plots/CITE_CD4-ridge.pdf"), 
                 width = 4, height = 5)


plot_me <- subset(coembed, tech4 =="PHAGE" & seurat_clusters != "7" & seurat_clusters != "13" & seurat_clusters != "12" )
DefaultAssay(plot_me) <- "PDT"
Idents(plot_me) <- plot_me$celltype
plot_me <- NormalizeData(plot_me, assay = "PDT", normalization.method = "CLR")

pgo2 <- RidgePlot(plot_me, features  = "NB17", sort = "increasing") + 
  ggtitle(NULL) +   scale_fill_manual(values = col_vec) +
  theme(legend.position = "none") 

cowplot::ggsave2(pgo2,
                 file = paste0("../plots/Phage_CD4-ridge.pdf"), 
                 width = 4, height = 5)

