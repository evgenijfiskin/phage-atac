library(Seurat)
library(data.table)
library(dplyr)
library(BuenColors)
library(Signac)

cells <- gsub("-1", "", fread("../data/barcodes.tsv", header = FALSE)[[1]])

ic <- function(donor_n){
  dt <- fread(paste0("../data/hashtag_",donor_n,"_barcodeCOUNTSfiltered.csv.gz"), header = FALSE)
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
dim(cmat)

counts <- Read10X_h5("../../../phage_atac_large_data_files/input/CD8_hashing/filtered_peak_bc_matrix.h5")
meta_data <- fread("../data/EF_hCD8_ATAC_singlecell.csv.gz") %>% filter(cell_id != "None")
barcode_vec <- meta_data[[1]]; meta_data <- data.frame(meta_data[,-1]); rownames(meta_data) <- barcode_vec
cd8s <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = meta_data
)

# Do deminsion reduction
cd8s <- RunTFIDF(cd8s)
cd8s <- FindTopFeatures(cd8s, min.cutoff = 'q0')
cd8s <- RunSVD(
  object = cd8s,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

# Make an embedding
cd8s <- RunUMAP(object = cd8s, reduction = 'lsi', dims = 2:20)
cd8s <- FindNeighbors(object = cd8s, reduction = 'lsi', dims = 2:20)
cd8s <- FindClusters(object = cd8s, resolution = 0.2, algorithm = 3)
DimPlot(object = cd8s, label = TRUE) + NoLegend()

cmat_for_seurat <- t(cmat); colnames(cmat_for_seurat) <- paste0(colnames(cmat_for_seurat), "-1")
cd8s[["ADT"]] <- CreateAssayObject(counts = cmat_for_seurat)
cd8s <- NormalizeData(cd8s, assay = "ADT", normalization.method = "CLR") %>% ScaleData(assay = "ADT")
cd8s$maxCD8 <- matrixStats::colMaxs(cd8s@assays$ADT@scale.data)
gene.activities <- readRDS("../../../phage_atac_large_data_files/output/CD8_hashing/CD8hashed.gene_activities.rds")
cd8s[["RNA"]] <- CreateAssayObject(counts = gene.activities)
FeaturePlot(cd8s, c("c50", "c51", "c54", "c55", "maxCD8"))
saveRDS(cd8s, file = "../../../phage_atac_large_data_files/output/CD8_hashing/CD8hashed.seuratDimreduction.rds")


