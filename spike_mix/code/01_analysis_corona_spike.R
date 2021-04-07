library(data.table)
library(Seurat)
library(Signac)
library(Matrix)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(viridis)
library(dplyr)
library(BuenColors)

# Don't rely on 10x knee call to see if we can recover a few more cells
raw_mat <- Read10X_h5("../data/raw_peak_bc_matrix.h5")
counts_mat <- raw_mat[,colSums(raw_mat) > 500]
dim(raw_mat); dim(counts_mat)

# Import single cell
metadata <- read.csv(
  file = "../data//singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts_mat,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '../../../phage_atac_large_data_files/input/spike-mix/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# Deal with annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
Annotation(pbmc) <- annotations

pbmc <- TSSEnrichment(pbmc)
pbmc$pct_mito <- pbmc$mitochondrial/pbmc$total *100

ggplot(pbmc@meta.data, aes(x = TSS.enrichment, y = pct_mito)) + 
  geom_point(size = 0.1)

pbmc <- subset(pbmc, TSS.enrichment >= 4)

# Skip additional filtering for now
# Jump straight to dimension reduction
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q50')
pbmc <- RunSVD(pbmc)

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3, resolution = 0.8)

# Import and process nanobody tags
dt2 <- fread("../data//EF_PBMC_293_PH.kallisto.counts.tsv.gz")
phagetags <- data.matrix(dt2[,-1])
rownames(phagetags) <- paste0(dt2[[1]], "-1")
ptfilt <- t(phagetags[colnames(pbmc),])
pbmc[["phage"]] <- CreateAssayObject(counts = ptfilt)
pbmc <- NormalizeData(pbmc, assay = "phage", normalization.method = "CLR")
pbmc <- ScaleData(pbmc, assay = "phage")

all_nbs <- c("NB17","NB18","NB19","NB25","NB44","NB65","NB66","NB67","NB72","NB75","NB85","NB86")

plotter <- c("NB25", "NB18", "NB17")
FeaturePlot(pbmc, features = plotter, 
            min.cutoff = "q05", max.cutoff = "q95", pt.size = 0.1,
            order = FALSE) &
  scale_color_viridis()
DimPlot(object = pbmc, label = TRUE) + NoLegend()

sc <- pbmc@meta.data$seurat_clusters
cd4 <- pbmc@assays$phage@data["NB17",]
cd8 <- pbmc@assays$phage@data["NB25",]
u1 <- pbmc@reductions$umap@cell.embeddings[,1]
u2 <- pbmc@reductions$umap@cell.embeddings[,2]

sc_cl <- case_when(
  u1 < 2 & u2 < -2 ~ "z293ts",
  sc == 11 ~ "DCs",
  sc == 3 ~ "Naive_Bcell",
  sc == 8 ~ "Memory_Bcell",
  sc %in% c(4,7) ~ "Monos",
  sc %in% c(0,9) ~ "z293ts",
  sc == 1 & cd8 > 1 ~ "Naive_CD8T",
  sc == 1 ~ "Naive_CD4T", 
  sc == 5 ~ "Effector_CD8T",
  sc == 6 & cd8 > 1 ~ "Effector_CD8T",
  sc == 2 & cd8 > 1 ~ "Effector_CD8T",
  sc %in% c(2,10) ~ "Mem_CD4T", 
  sc %in% c(6) ~ "NKcell", 
  TRUE ~ "other"
)

table(sc_cl)
pbmc$sc_cl <- sc_cl
DimPlot(object = pbmc, label = TRUE, group.by = "sc_cl", pt.size = 0.1) + NoLegend()

# Compute gene activities
gene.activities <- GeneActivity(pbmc)
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- 'peaks'

FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'NCAM1', 'LYZ', 'CD4', 'CD8B'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3, order = TRUE
)


saveRDS(pbmc, file = "../../../phage_atac_large_data_files/output/spike-mix/3March2021_Seurat_object.rds")
