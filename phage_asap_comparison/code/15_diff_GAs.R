library(Seurat)
library(BuenColors)
library(dplyr)

load("../../../phage_atac_large_data_files/output/asap_phage_comparison/15Feb2021_analyzed_seurat_objects.rda")
coembed@meta.data$barcode <- rownames(coembed@meta.data)
t_cells <- (coembed@meta.data %>% filter(tech4 == "PHAGE" & seurat_clusters %in% c("0", "1", "5", "6", "8")))%>% pull(barcode)
pbmc.atac <- subset(coembed, tech4 == "PHAGE")
pbmc.atac <- NormalizeData(pbmc.atac, assay = "PDT", normalization.method = "CLR")
mdf <- data.frame(
  barcode = colnames(pbmc.atac@assays$PDT@data),
  t(pbmc.atac@assays$PDT@data)
)

amat <- (mdf[,c(2:5)])
mdf$assign_CD4_CD8 <- case_when(
  !(mdf$barcode %in% t_cells) ~ "other",
  amat[,"NB17"] > 1.5 &  amat[,"NB25"] < 1.2 ~ "CD4",
  amat[,"NB17"] < 1.0 &  amat[,"NB25"] > 1.2 ~ "CD8",
  TRUE ~ "other")

table(mdf$assign_CD4_CD8 )
pbmc.atac@meta.data$assign_CD4_CD8 <- mdf$assign_CD4_CD8
pbmc.atac@meta.data$CD16assign <- mdf$CD16assign

tcelldf <- mdf %>% filter(barcode %in% t_cells)
tcelldf$density <- get_density(tcelldf$CD4, tcelldf$CD8)

ggplot(tcelldf , aes(x = CD4, y = CD8, color = assign_CD4_CD8)) + 
  geom_point(size = 0.5) + theme_classic() + 
  scale_color_manual(values = c("CD4"="#3373ba", "other" = "#c4c6c4", "CD8" = "#ec2027")) +
  theme(legend.position = "bottom") + labs(x = "CLR CD4", y = "CLR CD8", color = "assignment") 

  

# Now set idents
pbmc.atac <- SetIdent(pbmc.atac, value = "assign_CD4_CD8")
DefaultAssay(pbmc.atac) <- "ACTIVITY"
find_markers_T <- FindMarkers(pbmc.atac, ident.1 = "CD4", ident.2 = "CD8",
                  min.pct = 0.1, logfc.threshold = 0, print.bar = TRUE)
find_markers_T$gene <- rownames(find_markers_T)
write.table(find_markers_T, file = "../output/diff_GA_Tcell-26feb.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

