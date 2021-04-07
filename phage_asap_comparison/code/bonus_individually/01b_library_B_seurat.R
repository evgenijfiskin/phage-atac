library(Signac)
library(data.table)
library(Seurat)
library(Matrix)
library(dplyr)
library(BuenColors)
library(EnsDb.Hsapiens.v86)
library(viridis)

import_kite_counts <- function(path){
  mtx <- fread(paste0(path,"/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0(path,"/featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- (fread(paste0(path,"/featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}

# Import the protein count data
adtmat <- import_kite_counts("../data/LibB/")
colnames(adtmat) <- paste0(colnames(adtmat), "-1")
phage_df <- fread("../data/LibB/EF_PBMCB.kallisto.counts.tsv")
phagemat <- t(data.matrix(phage_df[,2:6]))
colnames(phagemat) <- paste0(phage_df[[1]], "-1")

# Import the ATAC data
counts <- Read10X_h5(filename = "../../EF_PBMC_B_ATAC_hg38_v12-mtMask/outs/filtered_peak_bc_matrix.h5")

# Slower but convenient for the row parsing
metadata <- read.csv(
  file = "../../EF_PBMC_B_ATAC_hg38_v12-mtMask/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

CA <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '../../EF_PBMC_B_ATAC_hg38_v12-mtMask/outs/fragments.tsv.gz',
  min.cells = 0,
  min.features = 0
)
Annotation(CA) <- annotations
pbmcs <- CreateSeuratObject(
  counts = CA,
  assay = "peaks",
  meta.data = metadata
)

# Add the protein data to the seurat object and normalize
pbmcs[["ADT"]] <-  CreateAssayObject(counts = adtmat[,colnames(pbmcs)])
pbmcs <- NormalizeData(pbmcs, assay = "ADT", normalization.method = "CLR")
pbmcs <- ScaleData(pbmcs, assay = "ADT")

pbmcs[["PDT"]] <-  CreateAssayObject(counts = phagemat[,colnames(pbmcs)])
pbmcs <- NormalizeData(pbmcs, assay = "PDT", normalization.method = "CLR")
pbmcs <- ScaleData(pbmcs, assay = "PDT")

# Do ATAC-seq dimension reduction
pbmcs <- RunTFIDF(pbmcs) %>% FindTopFeatures( min.cutoff = 'q50') %>% 
  RunSVD() %>%  RunUMAP( reduction = 'lsi', dims = 2:30) %>%
  FindNeighbors(reduction = 'lsi', dims = 2:30) %>%
  FindClusters( verbose = FALSE,  resolution = 0.6)

DimPlot(object = pbmcs, label = TRUE) + NoLegend()


FeaturePlot(pbmcs, features = rownames(pbmcs@assays$ADT), min.cutoff = "q05", max.cutoff = "q95") & scale_color_viridis()
  
FeaturePlot(pbmcs, features = rownames(pbmcs@assays$PDT), min.cutoff = "q05", max.cutoff = "q95") & scale_color_viridis()



