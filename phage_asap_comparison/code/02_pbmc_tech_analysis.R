library(data.table)
library(Seurat)
library(Signac)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(viridis)
library(Matrix)
#-----
# Analyze PBMC phage atac
#-----
peaks <- Read10X_h5("../data/aggr/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../data/aggr/singlecell.csv.gz",
  header = TRUE,
  row.names = 1
) %>% dplyr::filter(cell_id != "None")
head(metadata); dim(metadata)

chrom_assay <- CreateChromatinAssay(
  counts = peaks,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '../../../phage_atac_large_data_files/input/asap_comparison/fragments.tsv.gz',
  min.cells = 10,
  min.features = 10
)

pbmc.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(pbmc.atac) <- annotations

pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = readRDS("../../../phage_atac_large_data_files/output/asap_phage_comparison/asap_phage_comparison_aggr.gene.activities.rds"))
pbmc.atac$tech <- "atac"

DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- FindVariableFeatures(pbmc.atac)
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac)

DefaultAssay(pbmc.atac) <- "ATAC"
pbmc.atac <- RunTFIDF(pbmc.atac) %>% FindTopFeatures( min.cutoff = 'q50') %>% 
  RunSVD() %>%  RunUMAP( reduction = 'lsi', dims = 2:30) %>%
  FindNeighbors(reduction = 'lsi', dims = 2:30) %>%
  FindClusters( verbose = FALSE,  resolution = 0.6)

#-----
# Now process cite-seq data
#-----

# Import ADT counts
dt <- fread("../data/previous_pub/GSE100866_PBMC_vs_flow_10X-ADT_umi.csv.gz")
adt <- data.matrix(dt[,-1])
rownames(adt) <- dt[[1]]

# Filter RNA matrix for human reads / data
dt2 <- fread("../data/previous_pub/GSE100866_PBMC_vs_flow_10X-RNA_umi.csv.gz")
rna <- data.matrix(dt2[,-1])
rownames(rna) <- dt2[[1]]

stopifnot(all(colnames(rna) == colnames(adt)))

# Get rid of the mouse cells... easy to see with a quick visualization
mdf <- data.frame(
  h_total = colSums(rna[substr(rownames(rna), 1, 1) == "H",]),
  m_total = colSums(rna[substr(rownames(rna), 1, 1) == "M",])
)
rna_filt <- rna[substr(rownames(rna), 1, 1) == "H",mdf$m_total < 100]
rownames(rna_filt) <- gsub("HUMAN_", "", rownames(rna_filt))

adt_filt <- adt[,mdf$m_total < 100]

dim(rna_filt)
dim(rna_filt)

# Do Seurat things
pbmc.cite <- CreateSeuratObject(counts = rna_filt, assay = "RNA", project = "10x_RNA")
pbmc.cite <- NormalizeData(pbmc.cite)
pbmc.cite <- FindVariableFeatures(pbmc.cite)
pbmc.cite <- ScaleData(pbmc.cite)
pbmc.cite <- RunPCA(pbmc.cite, verbose = FALSE)
pbmc.cite <- FindNeighbors(pbmc.cite, dims = 1:25)
pbmc.cite <- FindClusters(pbmc.cite, resolution = 0.8)
pbmc.cite <- RunUMAP(pbmc.cite, dims = 1:25)
DimPlot(pbmc.cite, reduction = "umap")

# Now do label transfer / same embedding
transfer.anchors <- FindTransferAnchors(reference = pbmc.cite, query = pbmc.atac,
                                        features = VariableFeatures(object = pbmc.cite), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

# joint embedding
genes.use <- VariableFeatures(pbmc.cite)
refdata <- GetAssayData(pbmc.cite, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]], dims = 2:30)
pbmc.atac[["RNA"]] <- imputation
coembed <- merge(x = pbmc.cite, y = pbmc.atac)
coembed$tech <- ifelse(!is.na(coembed$tech), coembed$tech, "RNA")
coembed$tech4 <- case_when(
  coembed$tech  == "RNA" ~ "CITE", 
  substr(colnames(coembed), 18, 18) == "3" ~ "PHAGE",
  substr(colnames(coembed), 18, 18) == "2" ~ "ASAPandPHAGE",
  substr(colnames(coembed), 18, 18) == "1" ~ "ASAP"
)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)

coembed <- RunUMAP(coembed, dims = 1:25, seed.use = 2)
coembed <- FindNeighbors(coembed, dims = 1:25)
coembed <- FindClusters(coembed, resolution = 0.4)
DimPlot(coembed, group.by = 'seurat_clusters', split.by = "tech4",  label = TRUE)

DimPlot(coembed,   label = TRUE)
DimPlot(coembed, group.by = 'seurat_clusters', split.by = "tech",  label = TRUE)


library(Matrix)

import_phage_counts <- function(lib, out){
  dt <- fread(paste0("../data/",lib,"/EF_PBMC_",lib,".kallisto.counts.tsv"))
  dm <- data.matrix(dt[,c(-1,-3)])
  rownames(dm) <- paste0(dt[[1]], out)
  return(t(dm))
}

# Import ASAP data
import_kite_counts <- function(lib, out){
  mtx <- fread(paste0("../data/",lib,"/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0("../data/",lib,"/featurecounts.barcodes.txt"), header = FALSE)[[1]], out)
  colnames(matx) <- paste0("asap", fread(paste0("../data/",lib,"/featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}
A_ADT <- import_kite_counts("LibA", "-1")
B_ADT <- import_kite_counts("LibB", "-2")
B_PDT <- import_phage_counts("LibB", "-2")
C_PDT <- import_phage_counts("LibC", "-3")

PDT <- sapply(rownames(B_PDT), function(pdt_marker){
  vec <- c(B_PDT[pdt_marker, ], C_PDT[pdt_marker, ]); names(vec) <- c(colnames(B_PDT), colnames(C_PDT))
  vec[colnames(coembed)]
}) %>% t()
PDT[is.na(PDT)] <- 0; colnames(PDT) <- colnames(coembed)

ASAP <- sapply(rownames(B_ADT), function(adt_marker){
  vec <- c(B_ADT[adt_marker, ], A_ADT[adt_marker, ]); names(vec) <- c(colnames(B_ADT), colnames(A_ADT))
  vec[colnames(coembed)]
}) %>% t()
ASAP[is.na(ASAP)] <- 0; colnames(ASAP) <- colnames(coembed)


# Add the CITE-seq counts
ADTmat <- sapply(rownames(adt_filt), function(adt_marker){
  vec <- c(adt_filt[adt_marker, ]); names(vec) <- c(colnames(adt_filt))
  vec[colnames(coembed)]
}) %>% t()
ADTmat[is.na(ADTmat)] <- 0; colnames(ADTmat) <- colnames(coembed)

coembed[["ADT"]] <- CreateAssayObject(counts = ADTmat[c("CD4", "CD8", "CD16"),])
coembed <- NormalizeData(coembed, assay = "ADT", normalization.method = "CLR")
coembed <- ScaleData(coembed, assay = "ADT")

# Add the ASAP/phage counts
coembed[["ASAP"]] <- CreateAssayObject(counts = ASAP)
coembed <- NormalizeData(coembed, assay = "ASAP", normalization.method = "CLR")
coembed <- ScaleData(coembed, assay = "ASAP")

coembed[["PDT"]] <- CreateAssayObject(counts = PDT)
coembed <- NormalizeData(coembed, assay = "PDT", normalization.method = "CLR")
coembed <- ScaleData(coembed, assay = "PDT")

FeaturePlot(coembed, features = c( rownames(PDT)), 
            min.cutoff = "q05", max.cutoff = "q95", split.by = "tech4", pt.size = 0.1) & scale_color_viridis()

vec <- c(
  "0" = "Naive_CD4T",
  "1" = "Mem_CD4T",
  "2" = "CD14_Mono",
  "3" = "Naive_Bcell",
  "4" = "NK_cell",
  "5" = "Effector_CD8T",
  "6" = "Naive_CD8T",
  "7" = "CD14_Mono",
  "8" = "Naive_CD4T",
  "9" = "Memory_Bcell",
  "10" = "CD16_Mono",
  "11" = "DCs",
  "12" = "Progenitors",
  "13" = "Platlets"
)
coembed@meta.data$celltype <- vec[as.character(coembed@meta.data$seurat_clusters)]

# Scale each assay separately
pbmc.cite[["ADT"]] <- CreateAssayObject(counts = ADTmat[c("CD4", "CD8", "CD16"),colnames(pbmc.cite)])
pbmc.cite <- NormalizeData(pbmc.cite, assay = "ADT", normalization.method = "CLR")
pbmc.cite <- ScaleData(pbmc.cite, assay = "ADT")
pbmc.cite[["umapd"]] <- CreateDimReducObject(embeddings = coembed@reductions$umap@cell.embeddings[coembed$tech != "atac",], key = "UMAPd")
pbmc.cite$seurat_clusters <- coembed@meta.data$seurat_clusters[coembed$tech != "atac"]

FeaturePlot(pbmc.cite, features = c("adt_CD16", "adt_CD4", "adt_CD8"),
            min.cutoff = "q01", max.cutoff = "q99", ncol = 3, reduction = "umapd")
Idents(pbmc.cite) <- pbmc.cite$seurat_clusters

# Append stuf to ATAC object
pbmc.atac[["ASAP"]] <- CreateAssayObject(counts = ASAP[,colnames(pbmc.atac)])
pbmc.atac <- NormalizeData(pbmc.atac, assay = "ASAP", normalization.method = "CLR")
pbmc.atac <- ScaleData(pbmc.atac, assay = "ASAP")

pbmc.atac[["PDT"]] <- CreateAssayObject(counts = PDT[,colnames(pbmc.atac)])
pbmc.atac <- NormalizeData(pbmc.atac, assay = "PDT", normalization.method = "CLR")
pbmc.atac <- ScaleData(pbmc.atac, assay = "PDT")

pbmc.atac[["umapd"]] <- CreateDimReducObject(embeddings = coembed@reductions$umap@cell.embeddings[coembed$tech == "atac",], key = "UMAPd", assay = "peaks")
pbmc.atac$seurat_clusters <- coembed@meta.data$seurat_clusters[coembed$tech == "atac"]
Idents(pbmc.atac) <- pbmc.atac$seurat_clusters

save(coembed, pbmc.atac, pbmc.cite, file = "../../../phage_atac_large_data_files/output/asap_phage_comparison/15Feb2021_analyzed_seurat_objects.rda")

