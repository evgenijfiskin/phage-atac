library(data.table)
library(dplyr)
library(Matrix)

import_phage_counts <- function(lib, out){
  dt <- fread(paste0("../../phage_asap_comparison/data/",lib,"/EF_PBMC_",lib,".kallisto.counts.tsv"))
  dm <- data.matrix(dt[,c(-1,-3)])
  rownames(dm) <- paste0(dt[[1]], out)
  return((dm))
}

# Import ASAP data
import_kite_counts <- function(lib, out){
  mtx <- fread(paste0("../../phage_asap_comparison/data/",lib,"/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0("../../phage_asap_comparison/data/",lib,"/featurecounts.barcodes.txt"), header = FALSE)[[1]], out)
  colnames(matx) <- paste0("asap", fread(paste0("../../phage_asap_comparison/data/",lib,"/featurecounts.genes.txt"), header = FALSE)[[1]])
  return((matx))
}
adt <- rbind(import_kite_counts("LibA", "-1"),import_kite_counts("LibB", "-2"))
pdt <- rbind(import_phage_counts("LibB", "-2"),import_phage_counts("LibC", "-3"))
adtdf <- data.frame(barcode = rownames(adt), data.matrix(adt))
pdtdf <- data.frame(barcode = rownames(pdt), data.matrix(pdt))

sc <- fread("../../phage_asap_comparison/data/aggr/singlecell.csv.gz")
mdf <- merge(merge(sc, adtdf, by = "barcode",all = TRUE), pdtdf, by = "barcode", all = TRUE)
mdf <- mdf[!is.na(mdf$cell_id),]
mdf[is.na(mdf)] <- 0
write.table(mdf, file = "../out/ASAP_phage_comparison_singlecell_counts.csv", 
            sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
