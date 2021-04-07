library(data.table)
library(Seurat)
library(Signac)
library(Matrix)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(viridis)
library(dplyr)
library(BuenColors)

pbmc <- readRDS("../../../phage_atac_large_data_files/output/spike-mix/3March2021_Seurat_object.rds")
set.seed(10)
color_vec <- c("#003600", sample(jdb_palette("corona", 10))[2:10])
p1 <- DimPlot(pbmc, group.by = "sc_cl", pt.size = 0.1) +
  scale_color_manual(values = color_vec) +
  theme_void() + ggtitle("")
cowplot::ggsave2(p1, file = "../plots/8March2021_Mix_embedding.pdf", width = 4.5, height = 2.7)

make_feature_plot_nb <- function(feature_name){
  DefaultAssay(pbmc) <- "phage"
  p1 <- FeaturePlot(pbmc, features = (feature_name), cols = viridis(n=500), 
                    min.cutoff = "q05", max.cutoff = "q95", pt.size = 0.1,
                    order = FALSE) + theme_void() + theme(legend.position = "none")
  cowplot::ggsave2(p1,file = paste0("../plots/markers/", feature_name, "_NB.pdf"), 
                   width = 3, height = 3)
  feature_name
}

all_nbs <- c("NB17","NB18","NB19","NB25","NB44","NB65","NB66","NB67","NB72","NB75","NB85","NB86")
sapply(all_nbs, function(x){
  quantile(pbmc@assays$phage@data[x,],c(0.05, 0.95))
})

all_GA <-  c("KLRD1", "CD8A", "LEF1", "MS4A1", "MAFB", "IL3RA", "CD4", "CEACAM4", "FCGR3A")
sapply(all_GA, function(x){
  quantile(pbmc@assays$RNA@data[x,],c(0.01, 0.9))
})

make_feature_plot_ga <- function(feature_name){
  DefaultAssay(pbmc) <- "RNA"
  p2 <- FeaturePlot(pbmc, features = (feature_name), 
                    min.cutoff = "q01", max.cutoff = "q90", pt.size = 0.1, sort = FALSE,
                    order = FALSE) + theme_void() + theme(legend.position = "none")
  cowplot::ggsave2(p2, file = paste0("../plots/markers/", feature_name, "_GA.pdf"), 
                   width = 3, height = 3)
  feature_name
}

lapply(all_GA, make_feature_plot_ga)


# Geometric mean
df <- data.frame(t(pbmc@assays$phage@data), pbmc@meta.data)
df %>% 
  mutate(x = sqrt(NB19 * NB44), y = (NB65 * NB66 * NB67 * NB72 * NB75 * NB85 * NB86)^(1/7)) %>%
  dplyr::filter(sc_cl == "z293ts" & x > 1.2 & y < 1.5 & (x > y)) %>% dim()

df %>% 
  mutate(x = sqrt(NB19 * NB44), y = (NB65 * NB66 * NB67 * NB72 * NB75 * NB85 * NB86)^(1/7)) %>%
  dplyr::filter(sc_cl == "z293ts" & x < 1.2 & y > 1 & (x < y))%>% dim()

set.seed(2)
p1 <- ggplot(df, aes(x = sqrt(NB19 * NB44), y = (NB65 * NB66 * NB67 * NB72 * NB75 * NB85 * NB86)^(1/7), color = sc_cl)) + 
  geom_point(size = 0.4)  + scale_color_manual(values = color_vec)  + 
  pretty_plot() + L_border()  +  theme(legend.position = "none")  +
  labs(x = "Geometric mean GFP", y = "Geometric mean spike")
cowplot::ggsave2(p1, file = "../plots/8March2021_geometric_mea.pdf", width = 2.4, height = 2.4)

p11 <- ggplot(df %>% filter(sc_cl == "z293ts"), aes(x = sqrt(NB19 * NB44), y = (NB65 * NB66 * NB67 * NB72 * NB75 * NB85 * NB86)^(1/7), color = sc_cl)) + 
  geom_point(size = 0.4)  + scale_color_manual(values = jdb_palette("corona")[1])  + 
  pretty_plot() + L_border()  +  theme(legend.position = "none")  +
  labs(x = "Geometric mean GFP", y = "Geometric mean spike")
cowplot::ggsave2(p11, file = "../plots/8March2021_geometric_mean_only293ts.pdf", width = 2.4, height = 2.4)

df <- data.frame(t(pbmc@assays$phage@counts), pbmc@meta.data)
df %>% 
  mutate(x = sqrt(NB19 * NB44), y = (NB65 * NB66 * NB67 * NB72 * NB75 * NB85 * NB86)^(1/7)) %>%
  dplyr::filter(sc_cl == "z293ts" & x > 1.2 & y < 1.5 & (x > y)) %>% dim()

df %>% 
  mutate(x = sqrt(NB19 * NB44), y = (NB65 * NB66 * NB67 * NB72 * NB75 * NB85 * NB86)^(1/7)) %>%
  dplyr::filter(sc_cl == "z293ts" & x < 1.2 & y > 1 & (x < y))%>% dim()

set.seed(2)
p1 <- ggplot(shuf(df), aes(x = sqrt(NB19 * NB44), y = (NB65 * NB66 * NB67 * NB72 * NB75 * NB85 * NB86)^(1/7), color = sc_cl)) + 
  geom_point(size = 0.4)  + scale_color_manual(values = color_vec)  + 
  pretty_plot() + L_border()  +  theme(legend.position = "none")  +
  labs(x = "Geometric mean GFP", y = "Geometric mean spike")
cowplot::ggsave2(p1, file = "../plots/10March2021_geometric_mean_counts-s.pdf", width = 2.4, height = 2.4)

p11 <- ggplot(df %>% dplyr::filter(sc_cl == "z293ts"), aes(x = sqrt(NB19 * NB44), y = (NB65 * NB66 * NB67 * NB72 * NB75 * NB85 * NB86)^(1/7), color = sc_cl)) + 
  geom_point(size = 0.4)  + scale_color_manual(values = jdb_palette("corona")[1])  + 
  pretty_plot() + L_border()  +  theme(legend.position = "none")  +
  labs(x = "Geometric mean GFP", y = "Geometric mean spike")
cowplot::ggsave2(p11, file = "../plots/10March2021_geometric_mean_only293ts_counts.pdf", width = 2.4, height = 2.4)



library(dplyr)
pheatmap::pheatmap(round(cor(t(pbmc@assays$phage@scale.data)),2),
                   method = "pearson", color = viridis(500))

make_feature_scatter <- function(feature1, feature2){
  p1 <- FeatureScatter(pbmc,  pt.size = 0.3, group.by = "sc_cl",slot="counts",
                       feature1 =feature1, feature2 = feature2, col = color_vec) +
    theme(legend.position = "none") 
  cowplot::ggsave2(p1, 
                  filename = paste0("../plots/pointscatter/", feature1, "_", feature2, "_scatter.pdf"), width = 3, height = 3)
                
}
matcombs <- combn(all_nbs,2)
sapply(1:dim(matcombs)[2], function(i){
  make_feature_scatter(matcombs[1,i], matcombs[2,i])
})

dfpp <- data.frame(color = pbmc@meta.data$sc_cl, 
                   t(pbmc@assays$phage@counts))
make_feature_scatter_manual <- function(feature1, feature2){
  dfpp1 <- dfpp[,c(feature1, feature2, "color")]
  colnames(dfpp1) <- c("featurep1", "featurep2", "color")
  mmax <- max(c(dfpp1[[1]], dfpp1[[2]]))
  p1 <- ggplot(dfpp1, aes(x = featurep1, y = featurep2, color = color))+
    geom_point(size = 0.4) + scale_color_manual(values = color_vec) +
    scale_y_continuous(limits = c(0, mmax)) + scale_x_continuous(limits = c(0, mmax)) +
    pretty_plot(fontsize = 8) + L_border() + labs(x = feature1, y = feature2) +
    theme(legend.position = "none") 
  cowplot::ggsave2(p1, 
                   filename = paste0("../plots/pointscatter/", feature1, "_", feature2, "_scatter_manual.pdf"), width = 3, height = 3)
  
}
matcombs <- combn(all_nbs,2)
sapply(1:dim(matcombs)[2], function(i){
  make_feature_scatter_manual(matcombs[1,i], matcombs[2,i])
})


df2 <- data.frame(t(pbmc@assays$phage@counts), sc_cl = pbmc@meta.data$sc_cl)
df2 %>% dplyr::filter(sc_cl == "z293ts") %>%
  reshape2::melt(id.vars = "sc_cl") %>%
  dplyr::filter(variable == "z293ts") -> melt_df2

p1 <- ggplot(melt_df2 %>% dplyr::filter(!(variable %in% c("NB17", "NB18", "NB25"))), aes(x = variable, y = value)) + 
  geom_boxplot(outlier.shape = NA, width = 0.5) + scale_y_log10() + 
  pretty_plot(fontsize = 8) + L_border()  +  theme(legend.position = "none")  +
  labs(x = "nanobody", y = " Counts") 
p1

cowplot::ggsave2(p1, file = "../plots/boxplot_spike_gfp.pdf", 
                 width = 4.2, height = 2)

