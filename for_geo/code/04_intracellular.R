library(data.table)
library(dplyr)

import_phage_counts <- function(){
  dt <- fread(paste0("../../intracellular/data/PA_GFP_21Feb2021.kallisto.counts.tsv"))
  dm <- data.matrix(dt[,c(-1,-3)])
  rownames(dm) <- paste0(dt[[1]],"-1")
  return(t(dm))
}

pdt <- t(import_phage_counts())[,6]
sc <- fread("../../intracellular/data/singlecell.csv")
sc$EGFP_quant <- ifelse(is.na(pdt[sc[[1]]]), 0, pdt[sc[[1]]])
write.table(sc, file = "../out/Intracellular_singlecell_counts.csv", 
            sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
