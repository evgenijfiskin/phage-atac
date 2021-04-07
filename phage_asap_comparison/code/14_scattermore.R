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

make_feature_scatter <- function(feature_name1, feature_name2, assay, techuse, DS = FALSE){
  
  plot_me <- subset(coembed, tech4 == techuse)
  DefaultAssay(plot_me) <- assay
  Idents(plot_me) <- plot_me$celltype
  plot_me <- NormalizeData(plot_me, assay = assay, normalization.method = "CLR")
  
  pgo <- FeatureScatter(plot_me, feature1 = feature_name1, feature2 = feature_name2, pt.size = 0.3) + 
    ggtitle(NULL) +   scale_color_manual(values = col_vec) +
    theme(legend.position = "none") 
  
  cowplot::ggsave2(pgo,
                   file = paste0("../plots/feature_scatter/", feature_name1, "_", feature_name2, "_", techuse, "_featurescatter-count-ds.pdf"), 
                   width = 3, height = 3)
  feature_name1
}
make_feature_scatter("NB25", "NB18", "PDT", "PHAGE")
make_feature_scatter("NB25", "NB18", "PDT", "ASAPandPHAGE")


make_feature_scatter("CD4", "CD8", "ADT", "CITE")
make_feature_scatter("asapCD4", "asapCD8", "ASAP", "ASAPandPHAGE")
make_feature_scatter("asapCD4", "asapCD8", "ASAP", "ASAP")
make_feature_scatter("NB17", "NB25", "PDT", "ASAPandPHAGE")
make_feature_scatter("NB17", "NB25", "PDT", "PHAGE")


# Make a dedicated version for CITE-seq that removes the weird doublet cluster
plot_me <- subset(coembed, tech4 =="CITE" & seurat_clusters != "7" )
feature_name1 <- "CD4"; feature_name2 <- "CD8"; assay <- "ADT"
if(TRUE){
  set.seed(1)
  plot_me$ss <- 1:dim(plot_me)[2] %in% sample(1:dim(plot_me)[2], 1408)
  plot_me <- subset(plot_me, ss)
}
DefaultAssay(plot_me) <- assay
Idents(plot_me) <- plot_me$celltype
plot_me <- NormalizeData(plot_me, assay = assay, normalization.method = "CLR")

pgo <- FeatureScatter(plot_me, feature1 = feature_name1, feature2 = feature_name2, pt.size = 0.3) + 
  ggtitle(NULL) +   scale_color_manual(values = col_vec) +
  theme(legend.position = "none") 

cowplot::ggsave2(pgo,
                 file = paste0("../plots/feature_scatter/", feature_name1, "_", feature_name2, "_CITE_featurescatter-counts-ds.pdf"), 
                 width = 3, height = 3)

ddf <- data.frame(
  phage4= coembed@assays$PDT@data[c("NB17"),coembed$tech4 == "ASAPandPHAGE"],
  asap4= coembed@assays$ASAP@data[c("asap CD4"),coembed$tech4 == "ASAPandPHAGE"],
  ct  = coembed@meta.data$celltype[coembed$tech4 == "ASAPandPHAGE"]
)

pA <- ggplot(ddf %>% filter(ct %in% c("Naive_CD4T", "Mem_CD4T")), aes(x = asap4, y = phage4, color = ct)) +
  ggtitle(NULL) +   scale_color_manual(values = col_vec) +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") +
  geom_point(size = 0.3)  + scale_x_continuous(limits = c(0,2.7))
cowplot::ggsave2(pA, file = "../plots/costain_cd4_CLR_T.pdf", width = 2, height = 2)
