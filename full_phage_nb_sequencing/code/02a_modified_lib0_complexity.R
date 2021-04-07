library(data.table)
library(dplyr)
library(BuenColors)


estimateLibrarySize <- function(nTotal, nUnique){
  
  f <- function(x, c, n) {
    return(c / x - 1 + exp(-n / x))
  }
  
  m = 1
  M = 100
  
  nDuplicates <- (nTotal - nUnique) + 1 # charity to handle only unique reads observed
  
  if (nUnique > nTotal | (f(m * nUnique, nUnique, nTotal) < 0) | nUnique < 0 | nTotal < 0 | nDuplicates < 0) {
    message("Library size returns 0 -- invalid inputs; check this cell more closely")
    return(0)
  }
  
  while (f(M * nUnique, nUnique, nTotal) > 0) {
    M <- M*10.0
  }
  
  for(i in seq(0, 40)) {
    r <-  (m + M) / 2.0
    u <- f(r * nUnique, nUnique, nTotal);
    if (u == 0) {
      break
    } else if (u > 0) {
      m = r
    } else if (u < 0) {
      M = r
    }
  }
  
  return(round(nUnique * (m + M) / 2.0))
}


plib0 <- fread(paste0("../../../phage_atac_large_data_files/input/phage-evolution-sequencing/data/PLIB0_parsed.csv.gz"))
plib3 <- fread(paste0("../../../phage_atac_large_data_files/input/phage-evolution-sequencing/data/PLIB3_parsed.csv.gz"))

count_3 <- plib3[,.N,c("CDR1aa", "CDR2aa", "CDR3aa")]
count_3 <- count_3[count_3$N >=2, ]
lib3clones <- paste0(count_3[[1]], count_3[[2]], count_3[[3]])
lib0id <- paste0(plib0[["CDR1aa"]], plib0[["CDR2aa"]], plib0[["CDR3aa"]])

plib0f <- plib0[!(lib0id %in% lib3clones),]
dim(plib0f)
dim(plib0)

estimateLibrarySize(dim(plib0f)[1], (plib0f[,.N,c("CDR1aa", "CDR2aa", "CDR3aa")] %>% dim())[1])
estimateLibrarySize(dim(plib0)[1], (plib0[,.N,c("CDR1aa", "CDR2aa", "CDR3aa")] %>% dim())[1])
