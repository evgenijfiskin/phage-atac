library(Seurat)
library(BuenColors)
library(dplyr)
library(viridis)
load("../../../phage_atac_large_data_files/output/pbmcs/28July2020_analyzed_seurat_objects.rda")

# Import phage atac data
import_tag <- function(tag, cells){
  dt_F8 <- fread(paste0("../data/F8/",tag,"barcodeCOUNTSfiltered.csv.gz"))
  dt_F9 <- fread(paste0("../data/F9/",tag,"barcodeCOUNTSfiltered.csv.gz"))
  dt_F8$V3 <- paste0(dt_F8[[2]], "-1")
  dt_F9$V3 <- paste0(dt_F9[[2]], "-2")
  vec <- c(dt_F8[[1]], dt_F9[[1]]); names(vec) <- c(dt_F8[[3]], dt_F9[[3]])
  unname(vec[cells])
}

cite.pa <- data.matrix(data.frame(
  CD4 = import_tag("CD4", colnames(pbmc.atac)),
  CD8 = import_tag("CD8", colnames(pbmc.atac)),
  CD16 = import_tag("CD16", colnames(pbmc.atac)),
  GFP = import_tag("GFP", colnames(pbmc.atac))
)) %>% t()
cite.pa[is.na(cite.pa)] <- 0
colnames(cite.pa) <- colnames(pbmc.atac)
pbmc.atac[["ADT"]] <- CreateAssayObject(counts = cite.pa)
pbmc.atac <- NormalizeData(pbmc.atac, assay = "ADT", normalization.method = "CLR")
pbmc.atac <- ScaleData(pbmc.atac, assay = "ADT")

p1 <- FeaturePlot(pbmc.atac, "GFP", reduction = "umapd", pt.size = 0.1) +
  theme_void() + scale_color_viridis() + 
  ggtitle("") + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/GFP_UMAP.pdf",
                 width = 4, height = 4)

