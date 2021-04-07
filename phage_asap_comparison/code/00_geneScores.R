library(EnsDb.Hsapiens.v86)
library(Seurat)
library(Signac)
library(dplyr)

#-----
# Analyze PBMC phage atac
#-----
peaks <- Read10X_h5("../data/aggr/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../data/aggr/singlecell.csv.gz",
  header = TRUE,
  row.names = 1
) %>% filter(cell_id != "None")
head(metadata); dim(metadata)

chrom_assay <- CreateChromatinAssay(
  counts = peaks,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '../../../phage_atac_large_data_files/input/asap_comparison/fragments.tsv.gz',
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
