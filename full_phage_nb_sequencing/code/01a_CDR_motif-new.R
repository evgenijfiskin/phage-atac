library(data.table)
library(motifStack)
library(dplyr)
library(stringr)
library(xlsx)

aas <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

make_table_mat <- function(dense_mat){
  sapply(aas, function(aa){
    colSums(dense_mat == aa)
  }) %>% t() -> count_mat
  colnames(count_mat) <- 1:dim(count_mat)[2]
  count_mat
}

norm <- function(mat){
  round((t(t(mat)/colSums(mat)) * 100 ), 2)
}

process_lib_collapse <- function(lib_id){
  plib <- fread(paste0("../../../phage_atac_large_data_files/input/phage-evolution-sequencing/data/",lib_id,"_parsed.csv.gz"))
  sdf <- plib[,.N,c("CDR1nt", "CDR2nt", "CDR3nt", "CDR1aa", "CDR2aa", "CDR3aa", "CDR3class")]
  plib <- sdf
  
  # Make CDR1 plot
  cdr1_mat <- str_split_fixed(plib[["CDR1aa"]], "", 7) %>% make_table_mat
  cdr2_mat <- str_split_fixed(plib[["CDR2aa"]], "", 12) %>% make_table_mat
  cdr3s_mat <- str_split_fixed(plib[["CDR3aa"]][plib$CDR3class == "short"], "", 10) %>% make_table_mat
  cdr3m_mat <- str_split_fixed(plib[["CDR3aa"]][plib$CDR3class == "medium"], "", 14) %>% make_table_mat
  cdr3l_mat <- str_split_fixed(plib[["CDR3aa"]][plib$CDR3class == "long"], "", 18) %>% make_table_mat
  
  # Write percentages to xlsx
  fn <- paste0("../output/motifs/", lib_id, "_percent-new.xlsx")
  write.xlsx(data.frame(norm(cdr1_mat)), file=fn, sheetName="CDR1", row.names=TRUE)
  write.xlsx(data.frame(norm(cdr2_mat)), file=fn, sheetName="CDR2", append=TRUE, row.names=TRUE)
  write.xlsx(data.frame(norm(cdr3s_mat)), file=fn, sheetName="CDR3short", append=TRUE, row.names=TRUE)
  write.xlsx(data.frame(norm(cdr3m_mat)), file=fn, sheetName="CDR3medium", append=TRUE, row.names=TRUE)
  write.xlsx(data.frame(norm(cdr3l_mat)), file=fn, sheetName="CDR3long", append=TRUE, row.names=TRUE)
  
  # Make Motif plot
  pdf(paste0("../output/motifs/", lib_id, "_CDRmotifs.pdf"), width = 5, height = 3)
  new("pfm", mat=pcm2pfm(cdr1_mat), name= paste0("CDR1 ", lib_id),
      color=colorset(alphabet="AA",colorScheme="chemistry")) %>% plot()
  new("pfm", mat=pcm2pfm(cdr2_mat), name= paste0("CDR2 ", lib_id),
      color=colorset(alphabet="AA",colorScheme="chemistry")) %>% plot()
  new("pfm", mat=pcm2pfm(cdr3s_mat), name= paste0("CDR3 ", lib_id,  " short"),
      color=colorset(alphabet="AA",colorScheme="chemistry")) %>% plot()
  new("pfm", mat=pcm2pfm(cdr3m_mat), name= paste0("CDR3 ", lib_id, " medium"),
      color=colorset(alphabet="AA",colorScheme="chemistry")) %>% plot()
  new("pfm", mat=pcm2pfm(cdr3l_mat), name= paste0("CDR3 ", lib_id, " long"),
      color=colorset(alphabet="AA",colorScheme="chemistry")) %>% plot()
  dev.off()
}

process_lib_collapse("PLIB3")


process_lib_collapse("PLIB0")
process_lib_collapse("PLIB1")
process_lib_collapse("PLIB2")
process_lib_collapse("PLIB4")