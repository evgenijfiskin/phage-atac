library(data.table)
library(dplyr)

# Import single-cell data
tag_df <- rbind(
  fread("../data/evolved/Fig2_C5_barcodeCOUNTSfiltered_human_clean.csv.gz"),
  fread("../data/evolved/Fig2_C5_barcodeCOUNTSfiltered_mouse_clean.csv.gz"),
  fread("../data/evolved/Fig2_C5_barcodeCOUNTSfiltered_noncells_clean.csv.gz")
)
tag_df$barcode <- paste0(tag_df$V2, "-1")
commercial_df <- tag_df[,c(3,1)]; colnames(commercial_df) <- c("barcode", "eGFP_counts")

# Merge with single cell data
tenx_df <- fread("../data/evolved/Fig2_singlecell.csv.gz")
commercial_df_all <- left_join(tenx_df, commercial_df, by = "barcode")
commercial_df_all$eGFP_counts <- ifelse(!is.na(commercial_df_all$eGFP_counts), commercial_df_all$eGFP_counts, 0)

commercial_df_all$FRIP <- commercial_df_all$peak_region_fragments/commercial_df_all$passed_filters
ggplot(commercial_df_all %>% dplyr::filter(cell_id != "None"), aes(x = passed_filters_GRCh38, y = passed_filters_mm10, color = FRIP)) + 
  geom_point(size = 0.2) + scale_y_log10() + scale_x_log10() +
  pretty_plot(fontsize = 8) + L_border() + scale_color_viridis()


# Make Plots

phage$purity <- phage$peak_region_fragments_GRCh38 / (phage$peak_region_fragments_GRCh38 + phage$peak_region_fragments_mm10)
phage <- commercial_df_all %>% filter(is_GRCh38_cell_barcode == 1 & peak_region_fragments_GRCh38 > 5000 & FRIP > 0.2 & purity > 0.98) %>% mutate(classification = ifelse(eGFP_counts > 1000, "high", "low"))

pPhage2 <- ggplot(phage , aes(x = log10(eGFP_counts), fill = classification)) +
  geom_histogram(data = phage %>% filter(classification == "low"), aes(y = ..density..), alpha = 0.3) +
  geom_histogram(data = phage %>% filter(classification == "high"), aes(y = ..density..), alpha = 0.3) +
  geom_density(aes(color = classification), fill = NA) +
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("firebrick", "dodgerblue2")) +
  scale_fill_manual(values = c("firebrick", "dodgerblue2")) +
  ggtitle("Phage") + labs(x = "log10 eGFP (reads)") +
  theme(legend.position = "none") + theme(plot.title = element_text(size = 8))
pPhage2
qplot(log10(phage$eGFP_counts))

#cowplot::ggsave2(cowplot::plot_grid(pFacs, pPhage, ncol = 1), 
#                 file = "../plots/histograms_gfp2.pdf", width = 1.3, height = 2.2)


egfp_only <- commercial_df_all %>% filter(cell_id != "None") 
count_mat <- data.matrix(t(egfp_only$eGFP_counts)); colnames(count_mat) <- paste0("N", as.character(1:dim(count_mat)[2])); rownames(count_mat) <- "eGFP"
so <- count_mat %>% CreateSeuratObject(assay = "ADT") %>% NormalizeData( assay = "ADT", normalization.method = "CLR") %>% ScaleData
egfp_only$eGFP_CLR <- so@assays$ADT@scale.data[1,]
qplot(egfp_only%>% filter(is_GRCh38_cell_barcode == 1) %>% pull(eGFP_CLR)) +
  labs(x = "CLR - Human", y = "") + pretty_plot() + L_border()
