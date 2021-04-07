library(data.table)
library(Seurat)
library(Signac)
library(Matrix)
library(viridis)
library(dplyr)
library(BuenColors)

"%ni%" <- Negate("%in%")

raw <- readRDS("../../../phage_atac_large_data_files/output/spike-mix/3March2021_Seurat_object.rds")
k293s <- raw@meta.data$sc_cl == "z293ts"
pheatmap::pheatmap(round(cor(t(raw@assays$phage@scale.data[,k293s])),2),
                   method = "pearson", color = viridis(500))
