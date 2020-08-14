library(Seurat)
library(dplyr)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

cd8s <- readRDS("../../../phage_atac_large_data_files/output/CD8_hashing/CD8hashed.seuratDimreduction.rds")
frags <- fread("../../../phage_atac_large_data_files/input/CD8_hashing/fragments.tsv.gz")
cdf <- cd8s@meta.data
clusters <- as.character(unique(cdf$seurat_clusters))
cdf$cell_id <- rownames(cdf)
sapply(clusters, function(cluster){
  
  possible_ids <- cdf %>% filter(seurat_clusters == cluster) %>% pull(cell_id) %>% as.character()
  cluster_gr <- frags %>% filter(V4 %in% possible_ids) %>%
    setnames(c("chr", "start", "end", "V4", "PCRn")) %>%
    makeGRangesFromDataFrame()
  
  reads_coverage <- coverage(cluster_gr)/length(cluster_gr)*1000000
  export.bw(reads_coverage, con = paste0("../../../phage_atac_large_data_files/output/CD8_hashing/CD8hash-c", as.character(cluster), ".bw"))
  cluster
}) -> bulk2
