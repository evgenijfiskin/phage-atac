library(data.table)
library(Seurat)
library(Signac)
library(dplyr)
library(viridis)

#-----
# Analyze PBMC phage atac
#-----
peaks <- Read10X_h5("../data//filtered_peak_bc_matrix.h5")
pbmc.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")

meta <- read.table("../data//singlecell.csv",
                   sep = ",", header = TRUE, row.names = 1, 
                   stringsAsFactors = FALSE)
meta <- meta[colnames(pbmc.atac), ]
meta <- meta[meta$passed_filters> 1000,]
dim(meta)
summary(meta$passed_filters)
summary(meta$DNase_sensitive_region_fragments/meta$passed_filters)

pbmc.atac <- AddMetaData(pbmc.atac, metadata = meta)
pbmc.atac$tech <- "atac"

pbmc.atac <-  subset(pbmc.atac, subset = passed_filters >= 1000)
dim(pbmc.atac)
pbmc.atac <- RunTFIDF(pbmc.atac) %>% FindTopFeatures( min.cutoff = 'q50') %>% 
  RunSVD() %>%  RunUMAP( reduction = 'lsi', dims = 2:30) %>%
  FindNeighbors(reduction = 'lsi', dims = 2:30) %>%
  FindClusters( verbose = FALSE,  resolution = 0.6)

import_phage_counts <- function(){
  dt <- fread(paste0("PA_GFP_21Feb2021.kallisto.counts.tsv"))
  dm <- data.matrix(dt[,c(-1,-3)])
  rownames(dm) <- paste0(dt[[1]],"-1")
  return(t(dm))
}

pdt <- import_phage_counts()

PDT <- sapply(rownames(pdt), function(pdt_marker){
  vec <- c(pdt[pdt_marker, ]); names(vec) <- c(colnames(pdt))
  vec2 <- vec[colnames(pbmc.atac)]
  vec2[is.na(vec2)] <- 0
  vec2
}) %>% t()
colnames(PDT) <- colnames(pbmc.atac)

# Add the ADT counts
pbmc.atac[["PDT"]] <- CreateAssayObject(counts = PDT)
pbmc.atac <- NormalizeData(pbmc.atac, assay = "PDT", normalization.method = "CLR")
pbmc.atac <- ScaleData(pbmc.atac, assay = "PDT")

FeaturePlot(pbmc.atac, features = rownames(pbmc.atac@assays$PDT), order = TRUE,
            min.cutoff = "q05", max.cutoff = "q95") & scale_color_viridis()

FeatureScatter(pbmc.atac,feature1 = "EGFP-NB19", feature2 = "passed_filters" ) +
  scale_y_log10()

qplot(pbmc.atac@assays$PDT@scale.data["EGFP-NB19",], bins = 20) +
  qplot(log10(pbmc.atac@assays$PDT@counts["EGFP-NB19",]), bins = 20)

table(log10(pbmc.atac@assays$PDT@counts["EGFP-NB19",] )> 2.3)
      