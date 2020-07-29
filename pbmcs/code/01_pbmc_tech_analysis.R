library(data.table)
library(Seurat)
library(Signac)
library(dplyr)

#-----
# Analyze PBMC phage atac
#-----
peaks <- Read10X_h5("../../../phage_atac_large_data_files/input/pbmcs/phage_atac_pbmcs/outs/filtered_peak_bc_matrix.h5")
pbmc.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = readRDS("../../../phage_atac_large_data_files/output/pbmcs/signac_bothChannels_ga.rds"))
meta <- read.table("../../../phage_atac_large_data_files/input/pbmcs/phage_atac_pbmcs/outs/singlecell.csv",
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
VariableFeatures(pbmc.atac) <- names(which(Matrix::rowSums(pbmc.atac) > 100))
pbmc.atac <- RunLSI(pbmc.atac, n = 50, scale.max = NULL)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 1:50)
DimPlot(pbmc.atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")

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
pbmc.cite <- CreateSeuratObject(counts = rna_filt)
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

# transfer clusters
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.cite$seurat_clusters, 
                                     weight.reduction = pbmc.atac[["lsi"]])
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
hist(pbmc.atac$prediction.score.max)
abline(v = 0.5, col = "red")

# joint embedding
genes.use <- VariableFeatures(pbmc.cite)
refdata <- GetAssayData(pbmc.cite, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]])
pbmc.atac[["RNA"]] <- imputation
coembed <- merge(x = pbmc.cite, y = pbmc.atac)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

coembed$seurat_clusters_m <- ifelse(!is.na(coembed$seurat_clusters), coembed$seurat_clusters, coembed$predicted.id)
coembed$tech <- ifelse(!is.na(coembed$tech), coembed$tech, "RNA")

DimPlot(coembed, group.by = 'seurat_clusters_m', split.by = "tech")

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
  CD4 = import_tag("CD4", colnames(peaks)),
  CD8 = import_tag("CD8", colnames(peaks)),
  CD16 = import_tag("CD16", colnames(peaks))
)) %>% t()
cite.pa[is.na(cite.pa)] <- 0
colnames(cite.pa) <- colnames(peaks)
combined.adt <- cbind(adt_filt[c("CD4", "CD8", "CD16"),], 
                      (cite.pa))


# Add the ADT counts
coembed[["ADT"]] <- CreateAssayObject(counts = combined.adt)
coembed <- NormalizeData(coembed, assay = "ADT", normalization.method = "CLR")
coembed <- ScaleData(coembed, assay = "ADT")

# Scale each assay separately
pbmc.cite[["ADT"]] <- CreateAssayObject(counts = adt_filt[c("CD4", "CD8", "CD16"),])
pbmc.cite <- NormalizeData(pbmc.cite, assay = "ADT", normalization.method = "CLR")
pbmc.cite <- ScaleData(pbmc.cite, assay = "ADT")
pbmc.cite[["umapd"]] <- CreateDimReducObject(embeddings = coembed@reductions$umap@cell.embeddings[coembed$tech != "atac",], key = "UMAPd")
FeaturePlot(pbmc.cite, features = c("adt_CD16", "adt_CD4", "adt_CD8"),
            min.cutoff = "q01", max.cutoff = "q99", ncol = 3, reduction = "umapd")

pbmc.atac$predicted.id <- factor(as.character(pbmc.atac$predicted.id),levels = levels(pbmc.cite$seurat_clusters))
pbmc.atac[["ADT"]] <- CreateAssayObject(counts = cite.pa)
pbmc.atac <- NormalizeData(pbmc.atac, assay = "ADT", normalization.method = "CLR")
pbmc.atac <- ScaleData(pbmc.atac, assay = "ADT")
pbmc.atac[["umapd"]] <- CreateDimReducObject(embeddings = coembed@reductions$umap@cell.embeddings[coembed$tech == "atac",], key = "UMAPd", assay = "ATAC")

# Try it all together
coembed@assays$ADT@scale.data <- cbind(pbmc.cite@assays$ADT@scale.data, pbmc.atac@assays$ADT@scale.data )

FeaturePlot(coembed, features = c("adt_CD16", "adt_CD4", "adt_CD8"),
            min.cutoff = "q01", max.cutoff = "q99", ncol = 3)

save(coembed, pbmc.atac, pbmc.rna, file = "../../../phage_atac_large_data_files/output/pbmcs/28July2020_analyzed_seurat_objects.rda")

