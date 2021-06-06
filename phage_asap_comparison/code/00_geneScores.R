library(EnsDb.Hsapiens.v86)
library(Seurat)
library(Signac)
library(dplyr)

#-----
# Analyze PBMC phage atac
#-----
peaks <- Read10X_h5("../data/aggr/filtered_feature_bc_matrix.h5")[[2]]

chrom_assay <- CreateChromatinAssay(
  counts = peaks,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = '../data/aggr/atac_fragments.tsv.gz',
  min.cells = 10,
  min.features = 10
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(pbmc) <- annotations
gene.activities <- GeneActivity(pbmc)
saveRDS(gene.activities, "../../../phage_atac_large_data_files/output/asap_phage_comparison/asap_phage_comparison_aggr.gene.activities.rds")
