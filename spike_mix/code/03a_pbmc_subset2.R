library(data.table)
library(Seurat)
library(Signac)
library(Matrix)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(viridis)
library(dplyr)
library(BuenColors)

"%ni%" <- Negate("%in%")

raw <- readRDS("../../../phage_atac_large_data_files/output/spike-mix/3March2021_Seurat_object.rds")
pbmc <- subset(raw, seurat_clusters %ni% c(0,11))

df <- data.frame(CD4counts = pbmc@assays$phage@counts["NB17",], 
                 CD4CLR = pbmc@assays$phage@data["NB17",], 
                 CD8counts = pbmc@assays$phage@counts["NB25",], 
                 CD8CLR = pbmc@assays$phage@data["NB25",], 
                 pbmc@meta.data)
df$celltype3 <- ifelse(df$sc_cl %in% c("Mem_CD4T", "Naive_CD4T"), "CD4",
                       ifelse(df$sc_cl %in% c("Effector_CD8T", "Naive_CD8T"), "CD8",  "other"))

df[,c("celltype3", "CD4counts", "CD8counts")] %>% 
  filter( celltype3 %in% c("CD4", "CD8")) %>%
  mutate(what = ifelse(CD4counts > CD8counts, "CD4", "CD8")) %>%
  mutate(value = ifelse(CD4counts > CD8counts,CD4counts,CD8counts )) -> long_df
  
  
ggplot(long_df %>% filter(value > 10), aes(x = log10(value), fill = what)) +
  geom_density(alpha = 0.5)+ 
  scale_fill_manual(values = c("firebrick", "dodgerblue2", "forestgreen")) +
  pretty_plot() -> p1

cowplot::ggsave2(p1, file = "../plots/histo_dens.pdf", width = 3.5, height  = 3)
