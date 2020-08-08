library(data.table)
library(dplyr)

single_cell <- fread("../../cd8_hashing/data/EF_hCD8_ATAC_singlecell.csv.gz")
cells <- single_cell$barcode
ic <- function(donor_n){
  dt <- fread(paste0("../../cd8_hashing/data/hashtag_",donor_n,"_barcodeCOUNTSfiltered.csv.gz"))
  vec <- dt[[1]]; names(vec) <- paste0(dt[[2]], "-1")
  unname(vec[cells])
}

cmat <- data.matrix(data.frame(
  c50 = ic("50"),
  c51 = ic("51"),
  c54 = ic("54"),
  c55 = ic("55")
))
cmat[is.na(cmat)] <- 0

all_cd8 <- data.frame(
  single_cell,
  cmat
)

write.table(all_cd8, file = "../out/CD8_Human_T-cell_Hashing_singlecell_counts.csv", 
            sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
