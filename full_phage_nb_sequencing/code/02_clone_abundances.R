library(data.table)
library(dplyr)
library(BuenColors)
library(preseqR)


estimate_complexity_preseq <- function(clone_df){
  count_OG <- clone_df %>% group_by(N) %>% 
    summarize(count = n()) %>% mutate(count_norm = count/N) %>% arrange((N))
  count_OG_mat <- data.matrix(count_OG[,c("N", "count_norm")])
  original_generating_function <- preseqR.rSAC(count_OG_mat)
  return(original_generating_function(10^20))
}

count_clones_complexity_preseq <- function(lib_id, n_get = 1750000){
  print(lib_id)
  plib <- fread(paste0("../../data/",lib_id,"_parsed.csv.gz"))
  sapply(1:20, function(i){
    set.seed(i)
    idx <- sample(1:dim(plib)[1], n_get)
    estimate_complexity_preseq(plib[idx,.N,c("CDR1aa", "CDR2aa", "CDR3aa")])
  }) -> vec
  vec
}

# Function translated from java version: https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/sam/DuplicationMetrics.java
# Not vectorized!!! 
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

count_clones_complexity <- function(lib_id, n_get = 1750000){
  print(lib_id)
  plib <- fread(paste0("../../data/",lib_id,"_parsed.csv.gz"))
  sapply(1:20, function(i){
    set.seed(i)
    idx <- sample(1:dim(plib)[1], n_get)
    estimateLibrarySize(n_get, (plib[idx,.N,c("CDR1aa", "CDR2aa", "CDR3aa")] %>% dim())[1])
  }) -> vec
  vec
}

data.frame(
  P0 = count_clones_complexity_preseq("PLIB0"),
  P1 = count_clones_complexity_preseq("PLIB1"),
  P2 = count_clones_complexity_preseq("PLIB2"),
  P3 = count_clones_complexity_preseq("PLIB3")
) -> complexity_df

p1 <- ggplot(reshape2::melt(complexity_df), aes(x = variable, y = value)) + 
  geom_boxplot() + pretty_plot(fontsize = 8) + L_border() +
  scale_y_log10() + 
  labs(x = "Library", y = "Estimated Complexity")

cowplot::ggsave2(p1, file = "../output/guestimated_counts_preseq.pdf", width =2, height = 2)

