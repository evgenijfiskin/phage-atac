library(Seurat)
library(data.table)
library(dplyr)
library(BuenColors)
library(Signac)
library(EnsDb.Hsapiens.v86)
counts_mat <- Read10X_h5("../data/filtered_peak_bc_matrix.h5")

# Import single cell
metadata <- read.csv(
  file = "../data/EF_hCD8_ATAC_singlecell.csv.gz",
  header = TRUE,
  row.names = 1
) %>% dplyr::filter(cell_id != "None")
summary(metadata$DNase_sensitive_region_fragments/metadata$passed_filters)
summary(metadata$passed_filters)
dim(metadata)
chrom_assay <- CreateChromatinAssay(
  counts = counts_mat,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = "../../../phage_atac_large_data_files/input/CD8_hashing/fragments.tsv.gz", # dummy
  min.cells = 10,
  min.features = 200
)

cd8s <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# Deal with annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
Annotation(cd8s) <- annotations

cd8s <- TSSEnrichment(cd8s, fast = TRUE)
summary(cd8s$TSS.enrichment)
