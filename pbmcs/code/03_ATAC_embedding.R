library(Seurat)
library(Signac)

peaks <- Read10X_h5("../../../phage_atac_large_data_files/input/pbmcs/phage_atac_pbmcs_aggr/outs/filtered_peak_bc_matrix.h5")
pbmc.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = readRDS("../../../phage_atac_large_data_files/output/pbmcs/signac_bothChannels_ga.rds"))
meta <- read.table("../../../phage_atac_large_data_files/input/pbmcs/phage_atac_pbmcs_aggr/outs/singlecell.csv",
                   sep = ",", header = TRUE, row.names = 1, 
                   stringsAsFactors = FALSE)
meta <- meta[colnames(pbmc.atac), ]
pbmc.atac <- AddMetaData(pbmc.atac, metadata = meta)
pbmc.atac$tech <- "atac"

DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- FindVariableFeatures(pbmc.atac)
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac)

DefaultAssay(pbmc.atac) <- "ATAC"
VariableFeatures(pbmc.atac) <- names(which(Matrix::rowSums(pbmc.atac) > 0))
pbmc.atac <- RunLSI(pbmc.atac, n = 20, scale.max = NULL)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:20)
pbmc.atac <- FindNeighbors(pbmc.atac, reduction = "lsi")
pbmc.atac <- FindClusters(pbmc.atac, resolution = 0.2)
DimPlot(pbmc.atac, reduction = "umap", label = TRUE)  + ggtitle("scATAC-seq")
cluster9 <- rownames(pbmc.atac@meta.data)[pbmc.atac@meta.data$seurat_clusters == "9"]

FeaturePlot(object = pbmc.atac,
            features = c('CD3D','CD4','CD8A','LEF1','BCL11B','NKG7',
                         'KIR2DL4','KLRD1','ZEB2','IL3RA','HDC',
                         'CCR3','MAFB','TREM1','CD14','MS4A1','PAX5','CD38','FLT3'),
            pt.size = 0.1,max.cutoff ='q90',ncol = 5)
