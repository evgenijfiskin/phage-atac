library(data.table)

# Import phage atac data
import_tag <- function(tag, lib, cells){
  dt <- fread(paste0("../../pbmcs/data/",lib,"/",tag,"barcodeCOUNTSfiltered.csv.gz"))
  dt$V3 <- paste0(dt[[2]], "-1")
  vec <- c(dt[[1]]); names(vec) <- c(dt[[3]])
  unname(vec[cells])
}

c1_sc <- fread("../../pbmcs/data/F8/singlecell.csv")
c2_sc <- fread("../../pbmcs/data/F9/singlecell.csv")

c1_counts <- data.matrix(data.frame(
  CD4 = import_tag("CD4", "F8", c1_sc$barcode),
  CD8 = import_tag("CD8",  "F8", c1_sc$barcode),
  CD16 = import_tag("CD16",  "F8", c1_sc$barcode)
))

c1_counts[is.na(c1_counts)] <- 0

c2_counts <- data.matrix(data.frame(
  CD4 = import_tag("CD4", "F9", c2_sc$barcode),
  CD8 = import_tag("CD8",  "F9", c2_sc$barcode),
  CD16 = import_tag("CD16",  "F9", c2_sc$barcode)
))

c2_counts[is.na(c2_counts)] <- 0

c1_all <- data.frame(c1_sc, c1_counts)
c2_all <- data.frame(c2_sc, c2_counts)

write.table(c1_all, file = "../out/PBMC_Channel1_singlecell_counts.csv", 
            sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(c2_all, file = "../out/PBMC_Channel2_singlecell_counts.csv", 
            sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)


