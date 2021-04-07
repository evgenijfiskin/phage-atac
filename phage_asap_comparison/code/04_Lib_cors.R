library(Signac)
library(Seurat)
library(BuenColors)
library(Matrix)
library(pheatmap)
library(viridis)

load("../../../phage_atac_large_data_files/output/asap_phage_comparison/15Feb2021_analyzed_seurat_objects.rda")
libB <- subset(coembed, tech4 == "ASAPandPHAGE")

cormat <- cor(cbind(t(libB@assays$ASAP@scale.data),
                   t(libB@assays$PDT@scale.data)), method = "spearman")
round(cormat,3)
pheatmap(cormat, use = "pairwise.complete", color = viridis(500), fontsize = 6)

FeatureScatter(libB, "NB25", "asapCD8")
FeatureScatter(libB, "NB17", "asapCD4")
