library(Seurat)
library(data.table)
library(dplyr)
library(BuenColors)
library(Signac)

cd8s <- readRDS("../../../phage_atac_large_data_files/output/CD8_hashing/CD8hashed.seuratDimreduction_14Aug2020.rds")
p0 <- DimPlot(object = cd8s, label = FALSE, group.by = "seurat_clusters", size = 0.1) + 
  theme_void() + NoLegend() 
#cowplot::ggsave2(p0, file = "../plots/baseUmap.png", width = 4, height = 4, dpi = 500)

p1 <- FeaturePlot(cd8s, c("c50"), max.cutoff = "q90") +scale_color_viridis() +
  theme_void() + NoLegend() + ggtitle("")

p2 <- FeaturePlot(cd8s, c("c51"),  max.cutoff = "q90") +scale_color_viridis() +
  theme_void() + NoLegend()  + ggtitle("")

p3 <- FeaturePlot(cd8s, c("c54"),max.cutoff = "q90") + scale_color_viridis() +
  theme_void() + NoLegend()  + ggtitle("")

p4 <- FeaturePlot(cd8s, c("c55"), max.cutoff = "q90") + scale_color_viridis() +
  theme_void() + NoLegend()  + ggtitle("")

sapply(rownames(cd8s@assays$ADT@data), function(x){
  quantile(cd8s@assays$ADT@data[x,],c(0.00, 0.90, 0.95))
})



cowplot::ggsave2(p1, file = "../plots/adt_c50.png", width = 4, height = 4, dpi = 500)
cowplot::ggsave2(p2, file = "../plots/adt_c51.png", width = 4, height = 4, dpi = 500)
cowplot::ggsave2(p3, file = "../plots/adt_c54.png", width = 4, height = 4, dpi = 500)
cowplot::ggsave2(p4, file = "../plots/adt_c55.png", width = 4, height = 4, dpi = 500)

cowplot::ggsave2(p1, file = "../plots/adt_c50.pdf", width = 4, height = 4)
cowplot::ggsave2(p2, file = "../plots/adt_c51.pdf", width = 4, height = 4)
cowplot::ggsave2(p3, file = "../plots/adt_c54.pdf", width = 4, height = 4)
cowplot::ggsave2(p4, file = "../plots/adt_c55.pdf", width = 4, height = 4)



#######
table(cd8s$ADT_classification.global)
Idents(cd8s) <- "ADT_maxID"
cd8s@meta.data$total_ADT <- Matrix::colSums(cd8s@assays$ADT@counts)
ggplot(cd8s@meta.data %>% filter(passed_filters < 20000), aes(x = ADT_classification.global, y = total_ADT)) + 
  geom_boxplot() 


# Compare with genotypes
gc <- fread("../data/Gokcen_clusters.tsv.gz")
gc$barcodeM <- gsub("-1", "", gc$barcode)

cd8s@meta.data$barcode <- rownames(cd8s@meta.data)
mdf <- merge(gc, cd8s@meta.data, by.x = "barcode", by.y = "barcode")

mdf$assign2 <- case_when(grepl("/", mdf$assignment) ~ "Doublet", mdf$assignment == "3" ~ "c54", mdf$assignment == "0" ~ "c50",  mdf$assignment == "1" ~ "c51",
                         mdf$assignment == "2" ~ "c55", TRUE ~ "other")
mdf %>% group_by(assign2, seurat_clusters) %>% summarize(count = n()) %>%
  group_by(assign2) %>% 
  mutate(prop = count/sum(count)) -> df2

pA <- ggplot(df2 %>% filter(assign2 %in% c("Doublet","c50", "c51", "c54", "c55")), aes(x = assign2, y = prop, fill = seurat_clusters)) +
  geom_bar(stat = "identity") + pretty_plot(fontsize = 8) + L_border() + labs(x = "Genotype Assignment", y = "Proportion") +
  theme(legend.position = "none")

cd8s@meta.data %>% group_by(hash.ID, seurat_clusters) %>% summarize(count = n()) %>%
  group_by(hash.ID) %>% 
  mutate(prop = count/sum(count)) -> df
pB <- ggplot(df, aes(x = hash.ID, y = prop, fill = seurat_clusters)) + 
  geom_bar(stat = "identity") + pretty_plot(fontsize = 8) + L_border() + labs(x = "Hash Assignment", y = "Proportion", fill = "") +
  theme(legend.position = "none") 

cowplot::ggsave2(cowplot::plot_grid(pA, pB, nrow = 1, rel_widths = c(0.9, 1.1)),
                 file = "../plots/proportions_slack_go.pdf", width = 4, height = 2)

FeaturePlot(cd8s, c("LEF1", "IKZF2", "GZMK", "CCL5"), min.cutoff = "q0", max.cutoff = "q90")
