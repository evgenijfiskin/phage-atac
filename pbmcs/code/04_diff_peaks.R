library(Seurat)
library(BuenColors)
library(dplyr)

load("../../../phage_atac_large_data_files/output/pbmcs/28July2020_analyzed_seurat_objects.rda")

ggplot(data.frame(t(pbmc.atac@assays$ADT@scale.data)), aes(x = CD4, y = CD8)) + 
  geom_point(size = 0.5)

amat <- t(pbmc.atac@assays$ADT@scale.data)
pbmc.atac@meta.data$assign_CD4 <- case_when(
  amat[,"CD4"] > 0 &  amat[,"CD8"] < 0.5 ~ "CD4",
  amat[,"CD4"] < 0 &  amat[,"CD8"] > 0.5 ~ "CD8",
  TRUE ~ "other")

ddf <- data.frame(t(pbmc.atac@assays$ADT@scale.data), assign_CD4 = pbmc.atac@meta.data$assign_CD4)
ggplot(ddf, aes(x = CD4, y = CD8, color = assign_CD4)) + 
  geom_point(size = 0.5)

# Now set idents
pbmc.atac <- SetIdent(pbmc.atac, value = "assign_CD4")
hm <- FindMarkers(pbmc.atac, ident.1 = "CD4", ident.2 = "CD8", min.pct = 0.1, logfc.threshold = 0, print.bar = TRUE)

