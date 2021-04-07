library(Signac)
library(Seurat)
library(BuenColors)
library(Matrix)
library(pheatmap)
library(viridis)
library(dplyr)
# Import data and subset to interesting variants from previous scripts
load("../../../phage_atac_large_data_files/output/asap_phage_comparison/15Feb2021_analyzed_seurat_objects.rda")
af_plot <-  readRDS("../output/Phage_ASAP_visualize_coverage.rds")
interesting_df <- readRDS("../output/Phage_ASAP_everything_meta.rds")
interesting_df <- interesting_df %>% filter(variant %in% rownames(af_plot))
vars_go <- interesting_df %>% arrange(desc(mean) ) %>% head(10) %>% pull(variant) 
hq_cells <- intersect(colnames(coembed), colnames(af_plot))

df <- data.frame(
  coembed@meta.data[hq_cells,],
  coembed@reductions$umap@cell.embeddings[hq_cells,],
  t(data.matrix(af_plot[vars_go,hq_cells])))

vars <- colnames(df)[3:12]
lapply(vars, function(i){
  df2 <- df[,c("UMAP_1", "UMAP_2")]
  df2$var <- df[,as.character(i)] > 0.1
  pMain1 <- ggplot(df2 %>% arrange((var)), aes(x = UMAP_1, y = UMAP_2, color = var )) + 
    geom_point(size = 0.2)  + theme_void() + ggtitle(i) + theme(legend.position = "none") + scale_color_manual(values = c("lightgrey", "blue")) 
  
  cowplot::ggsave2(pMain1, 
                   file = paste0("../plots/mito_muts/mito_mut_umap_10percentCutoff_", i, ".pdf"), height = 2, width =1.8)
  
})

pMain1 <- ggplot(df %>% dplyr::filter(tech4 == "ASAPandPHAGE") %>%  arrange((X16243T.C)), aes(x = UMAP_1, y = UMAP_2, color = X16243T.C *100)) + 
  geom_point(size = 0.2)  + theme_void() + ggtitle("X16243T.C %") + theme(legend.position = "bottom") + 
  scale_color_gradientn(colours = c("lightgrey", "firebrick")) + labs(color = "16243T>C")

cowplot::ggsave2(pMain1, 
                 file = paste0("../plots/mito_muts/mito_mut_umap_X152TCcolor-onlyCombinedChannel.pdf"), height = 3, width =2.2)

pMain2 <- ggplot(df %>% dplyr::filter(tech4 == "ASAPandPHAGE") %>%  arrange((X152T.C)), aes(x = UMAP_1, y = UMAP_2, color = X152T.C *100)) + 
  geom_point(size = 0.2)  + theme_void() + ggtitle("X152T.C %") + theme(legend.position = "bottom") + 
  scale_color_gradientn(colours = c("lightgrey", "firebrick"),  limits = c(0,10), oob = scales::squish) + labs(color = "152T>C")

cowplot::ggsave2(pMain2, 
                 file = paste0("../plots/mito_muts/mito_mut_umap_X152TCcolor-onlyCombinedChannel.pdf"), height = 3, width =2.2)

pMain3 <- ggplot(df %>% dplyr::filter(tech4 == "ASAPandPHAGE") %>%  arrange((X195T.C)), aes(x = UMAP_1, y = UMAP_2, color = X195T.C *100)) + 
  geom_point(size = 0.2)  + theme_void() + ggtitle("X195T.C %") + theme(legend.position = "bottom") + 
  scale_color_gradientn(colours = c("lightgrey", "firebrick"), limits = c(0,10), oob = scales::squish) + labs(color = "195T>C")

cowplot::ggsave2(pMain3, 
                 file = paste0("../plots/mito_muts/mito_mut_umap_X195TCcolor-onlyCombinedChannel.pdf"), height = 3, width =2.2)

pMain4 <- ggplot(df %>% dplyr::filter(tech4 == "ASAPandPHAGE") %>%  arrange((X3244G.A)), aes(x = UMAP_1, y = UMAP_2, color = X3244G.A *100)) + 
  geom_point(size = 0.2)  + theme_void() + ggtitle("X3244G.A %") + theme(legend.position = "bottom") + 
  scale_color_gradientn(colours = c("lightgrey", "firebrick"), limits = c(0,10), oob = scales::squish) + labs(color = "3244G>A")

cowplot::ggsave2(pMain4, 
                 file = paste0("../plots/mito_muts/mito_mut_umap_X3244GAcolor-onlyCombinedChannel.pdf"), height = 3, width =2.2)
