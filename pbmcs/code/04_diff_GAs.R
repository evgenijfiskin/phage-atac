library(Seurat)
library(BuenColors)
library(dplyr)

load("../../../phage_atac_large_data_files/output/pbmcs/28July2020_analyzed_seurat_objects.rda")
coembed@meta.data$barcode <- rownames(coembed@meta.data)
t_cells <- (coembed@meta.data %>% filter(orig.ident == "10x_ATAC" & seurat_clusters %in% c("1", "2", "5", "7", "8")))%>% pull(barcode)
monocytes <- (coembed@meta.data %>% filter(orig.ident == "10x_ATAC" & seurat_clusters %in% c("0", "6"))) %>% pull(barcode)
mdf <- data.frame(
  barcode = colnames(pbmc.atac@assays$ADT@scale.data),
  t(pbmc.atac@assays$ADT@scale.data)
)

amat <- (mdf[,c(2:4)])
mdf$assign_CD4_CD8 <- case_when(
  !(mdf$barcode %in% t_cells) ~ "other",
  amat[,"CD4"] > 0 &  amat[,"CD8"] < 0.5 ~ "CD4",
  amat[,"CD4"] < 0 &  amat[,"CD8"] > 0.3 ~ "CD8",
  TRUE ~ "other")
mdf$CD16assign <- case_when(
  !(mdf$barcode %in% monocytes) ~ "other",
  amat[,"CD16"] > 0.5 ~ "CD16high",
  TRUE ~ "CD16low")
table(mdf$assign_CD4_CD8 )
pbmc.atac@meta.data$assign_CD4_CD8 <- mdf$assign_CD4_CD8
pbmc.atac@meta.data$CD16assign <- mdf$CD16assign

tcelldf <- mdf %>% filter(barcode %in% t_cells)
tcelldf$density <- get_density(tcelldf$CD4, tcelldf$CD8)
#tcelldf$density <- ifelse(tcelldf$density > 0.2, 0.2, tcelldf$density)
ggplot(tcelldf  %>% filter(CD4 < 3 & CD8 < 4), aes(x = CD4, y = CD8, color = density)) + 
  geom_point(size = 0.5) + theme_classic() + 
  scale_color_viridis() + 
  scale_x_continuous(limits = c(-2.5,3)) +   scale_y_continuous(limits = c(-2.5,4)) +
  theme(legend.position = "bottom") + labs(x = "CLR CD4", y = "CLR CD8", color = "denisty") 

ggplot(tcelldf , aes(x = CD4, y = CD8, color = assign_CD4_CD8)) + 
  geom_point(size = 0.5) + theme_classic() + 
  scale_color_manual(values = c("CD4"="#3373ba", "other" = "#c4c6c4", "CD8" = "#ec2027")) +
  theme(legend.position = "bottom") + labs(x = "CLR CD4", y = "CLR CD8", color = "assignment") 

  

# Now set idents
pbmc.atac <- SetIdent(pbmc.atac, value = "assign_CD4_CD8")
DefaultAssay(pbmc.atac) <- "ACTIVITY"
#find_markers_T <- FindMarkers(pbmc.atac, ident.1 = "CD4", ident.2 = "CD8",
#                  min.pct = 0.1, logfc.threshold = 0, print.bar = TRUE)
#find_markers_T$gene <- rownames(find_markers_T)
# write.table(find_markers_T, file = "../output/diff_GA_Tcell.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Make monocyte density plot
mono_df <- mdf %>% filter(barcode %in% monocytes)
mono_df$CD16assign <- mono_df$CD16 > 0.5
mono_dens <- ggplot(mono_df, aes(x = CD16, fill = CD16assign)) + 
  geom_histogram(data = mono_df %>% filter(CD16assign == FALSE), aes(y = ..density..), bins = 50) +
  geom_histogram(data = mono_df %>% filter(CD16assign == TRUE), aes(y = ..density..), bins = 50) +
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("#ec6f2f", "#f49c61")) +
  scale_fill_manual(values = c("#ec6f2f", "#f49c61")) +
  labs(x = "CLR CD16") + geom_vline(xintercept = 0.5, linetype = 2) +
  theme(legend.position = "none") + theme(plot.title = element_text(size = 8))
cowplot::ggsave2(mono_dens, 
                 file = "../plots/mono_histogram.pdf", width = 1.5, height = 1.5)

pbmc.atac <- SetIdent(pbmc.atac, value = "CD16assign")
DefaultAssay(pbmc.atac) <- "ACTIVITY"
table(Idents(pbmc.atac))
#find_markers_mono <- FindMarkers(pbmc.atac, ident.1 = "CD16high", ident.2 = "CD16low",
#                              min.pct = 0.1, logfc.threshold = 0, print.bar = TRUE)
qdf <- data.frame(
  FCGR3A = pbmc.atac@assays$ACTIVITY@counts["FCGR3A",],
  pbmc.atac@meta.data)

qdf %>% group_by(CD16assign) %>%
  summarize(mean(FCGR3A))


