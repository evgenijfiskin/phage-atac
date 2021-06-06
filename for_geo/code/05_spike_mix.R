library(data.table)
library(dplyr)

import_phage_counts <- function(){
  dt <- fread(paste0("../../spike_mix/data/EF_PBMC_293_PH.kallisto.counts.tsv.gz"))
  dm <- data.matrix(dt[,c(-1)])
  rownames(dm) <- paste0(dt[[1]],"-1")
  return(t(dm))
}

pdt <- t(import_phage_counts())
pdt <- data.frame(pdt, barcode = rownames(pdt))
metadata <- fread(file = "../../spike_mix/data/singlecell.csv")
mdf <- merge(metadata, data.frame(pdt), by.x = "barcode", by.y = "barcode", all = TRUE)
mdf[is.na(mdf)] <- 0
write.table(mdf, file = "../out/SpikeMix_singlecell_counts.csv", 
            sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
