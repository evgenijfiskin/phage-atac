library(Seurat)
library(data.table)
library(dplyr)
library(BuenColors)

cells <- gsub("-1", "", fread("../data/barcodes.tsv")[[1]])

ic <- function(donor_n){
  dt <- fread(paste0("../data/hashtag_",donor_n,"_barcodeCOUNTSfiltered.csv.gz"))
  vec <- dt[[1]]; names(vec) <- dt[[2]]
  unname(vec[cells])
}

count_df <- data.frame(
  c50 = ic("50"),
  c51 = ic("51"),
  c54 = ic("54"),
  c55 = ic("55")
)
rownames(count_df) <- cells
cmat <- data.matrix(count_df)
cmat[is.na(cmat)] <- 0

# do seurat things
pbmc.hashtag <- CreateSeuratObject(counts = t(cmat), assay = "HTO")
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.995)
table(pbmc.hashtag$HTO_classification.global)
Idents(pbmc.hashtag) <- "HTO_maxID"
pR <- RidgePlot(pbmc.hashtag, assay = "HTO", 
                features = rownames(pbmc.hashtag[["HTO"]])[1:4], ncol = 4, cols = jdb_palette("corona"))
ggsave(pR, file = "../plots/HTO_ridges.pdf", width = 10, height = 2.5)  
FeatureScatter(pbmc.hashtag, feature1 = "c50", feature2 = "c51")

# Visualize
Idents(pbmc.hashtag) <- "HTO_classification.global"
pbmc.hashtag.subset <- subset(pbmc.hashtag, idents = "Negative", invert = TRUE)

# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = pbmc.hashtag.subset, assay = "HTO"))))
# Calculate tSNE embeddings with a distance matrix
pbmc.hashtag.subset <- RunTSNE(pbmc.hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 50)
pBase <- DimPlot(pbmc.hashtag.subset, pt.size = 0.01) +
  labs(x = "t-SNE1", y = "t-SNE2") +
  scale_color_manual(values = c("red2", "black")) + 
  theme_void() + theme(legend.position = "bottom") + ggtitle("")
ggsave(pBase, file = "../plots/doublets_tSNE.png", dpi = 500, width = 2, height = 2.3)

p1 <- FeaturePlot(pbmc.hashtag.subset, "c50", max.cutoff = 'q99', pt.size = 0.01) +
  labs(x = "t-SNE1", y = "t-SNE2") + theme_void() + theme(legend.position = "none") + ggtitle("")
p2 <- FeaturePlot(pbmc.hashtag.subset, "c51", max.cutoff = 'q99', pt.size = 0.01) +
  labs(x = "t-SNE1", y = "t-SNE2") + theme_void() + theme(legend.position = "none") + ggtitle("")
p3 <- FeaturePlot(pbmc.hashtag.subset, "c54", max.cutoff = 'q99', pt.size = 0.01) +
  labs(x = "t-SNE1", y = "t-SNE2") + theme_void() + theme(legend.position = "none") + ggtitle("")
p4 <- FeaturePlot(pbmc.hashtag.subset, "c55", max.cutoff = 'q99', pt.size = 0.01) +
  labs(x = "t-SNE1", y = "t-SNE2") + theme_void() + theme(legend.position = "none") + ggtitle("")
ggsave(cowplot::plot_grid(p1, p2, p3, p4, nrow = 2), file = "../plots/HTO_tSNE.png", dpi = 500, width = 4, height = 4)

# Make the heatmap
p1 <- HTOHeatmap(pbmc.hashtag, assay = "HTO", ncells = 5000) +
  theme(legend.position = "bottom") + scale_fill_gradientn(colors = c("purple",  "black", "yellow"))
ggsave(p1, file = "../plots/visualize_HTO_heatmap.png", width = 3, height = 2.3)


# Compare with Gokcen
gc <- fread("../data/Gokcen_clusters.tsv.gz")
gc$barcodeM <- gsub("-1", "", gc$barcode)
hto_df <- pbmc.hashtag@meta.data
hto_df$barcodeM <- rownames(hto_df)
mdf <- merge(gc, hto_df)

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
