library(Seurat)
library(data.table)
library(dplyr)
library(BuenColors)
library(viridis)

cd8s <- readRDS(file = "../../../phage_atac_large_data_files/output/CD8_hashing/CD8hashed.seuratDimreduction_14Aug2020.rds")
Idents(cd8s) <- "ADT_maxID"
pR <- RidgePlot(cd8s, assay = "ADT", 
                features = rownames(cd8s[["ADT"]])[1:4], ncol = 4, cols = jdb_palette("corona"))
ggsave(pR, file = "../plots/HTO_ridges.pdf", width = 10, height = 2.5)  

# Visualize
Idents(cd8s) <- "ADT_classification.global"
cd8s.hashtag.subset <- subset(cd8s, idents = "Negative", invert = TRUE)

# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = cd8s.hashtag.subset, assay = "ADT"))))
# Calculate tSNE embeddings with a distance matrix
cd8s.hashtag.subset <- RunTSNE(cd8s.hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 75)
pBase <- DimPlot(cd8s.hashtag.subset, pt.size = 0.01, reduction = "tsne") +
  labs(x = "t-SNE1", y = "t-SNE2") +
  scale_color_manual(values = c("red2", "dodgerblue3")) + 
  theme_void() + theme(legend.position = "bottom") + ggtitle("")
pBase
cowplot::ggsave2(pBase, file = "../plots/doublets_tSNE.pdf", width = 2, height = 2.3)

p1 <- FeaturePlot(cd8s.hashtag.subset, "c50", max.cutoff = 'q99', pt.size = 0.01, reduction = "tsne") +
  labs(x = "t-SNE1", y = "t-SNE2") + theme_void() + theme(legend.position = "none") + ggtitle("") + scale_color_viridis()
p2 <- FeaturePlot(cd8s.hashtag.subset, "c51", max.cutoff = 'q99', pt.size = 0.01, reduction = "tsne") +
  labs(x = "t-SNE1", y = "t-SNE2") + theme_void() + theme(legend.position = "none") + ggtitle("") + scale_color_viridis()
p3 <- FeaturePlot(cd8s.hashtag.subset, "c54", max.cutoff = 'q99', pt.size = 0.01, reduction = "tsne") +
  labs(x = "t-SNE1", y = "t-SNE2") + theme_void() + theme(legend.position = "none") + ggtitle("") + scale_color_viridis()
p4 <- FeaturePlot(cd8s.hashtag.subset, "c55", max.cutoff = 'q99', pt.size = 0.01, reduction = "tsne") +
  labs(x = "t-SNE1", y = "t-SNE2") + theme_void() + theme(legend.position = "none") + ggtitle("") + scale_color_viridis()
cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, p4, nrow = 1), file = "../plots/HTO_tSNE.pdf",  width = 8, height = 2)

# Make the heatmap
p1 <- HTOHeatmap(cd8s, assay = "ADT", ncells = 5000) +
  theme(legend.position = "bottom") + scale_fill_gradientn(colors = c("purple",  "black", "yellow"))
cowplot::ggsave2(p1, file = "../plots/visualize_HTO_heatmap.pdf", width = 3, height = 2.3)
cowplot::ggsave2(p1, file = "../plots/visualize_HTO_heatmap.png", width = 3, height = 2.3)


# Compare with Gokcen
gc <- fread("../data/Gokcen_clusters.tsv.gz")
gc$barcodeM <- gsub("-1", "", gc$barcode)
hto_df <- cd8s@meta.data
hto_df$barcode <- rownames(hto_df)
mdf <- merge(gc, hto_df, by = "barcode")

mdf %>% group_by(hash.ID, assignment) %>%
  summarize(count = n()) -> confusion_df

ca <- confusion_df$assignment
confusion_df$assign_genotype <- case_when(
  ca == "0" ~ "c50",
  ca == "1" ~ "c51",
  ca == "3" ~ "c54",
  ca == "2" ~ "c55",
  TRUE ~ "Doublet"
)
confusion_df2 <- confusion_df %>% filter(hash.ID != "Negative") %>%
  group_by(hash.ID, assign_genotype) %>% summarize(total = sum(count))
levels <- c("c50", "c51", "c54", "c55", "Doublet")
confusion_df2$hash.ID <- factor(as.character(confusion_df2$hash.ID ), levels)
confusion_df2$assign_genotype <- factor(as.character(confusion_df2$assign_genotype ), rev(levels))
confusion_df2 <- confusion_df2 %>% group_by(assign_genotype) %>% mutate(fill_val = total / sum(total))

p1 <- ggplot(confusion_df2, aes(x = hash.ID, y = assign_genotype, label = total, fill = fill_val*100, color = total > 1000)) +
  geom_tile(color = "black") +
  geom_text(size = 2.5) + pretty_plot(fontsize = 8) + scale_fill_gradientn(colors = jdb_palette("brewer_heat")) +
  theme(legend.position = "bottom") + 
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  scale_color_manual(guide = "none", values = c("black", "white")) + 
  labs(x = "Hash ID", y = "Cell genotype", fill = "% of genotypes")
cowplot::ggsave2(p1, file = "../plots/heatmap_hash_concordance.pdf", width = 2, height = 2.5)
