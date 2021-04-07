library(BuenColors)
library(dplyr)
library(data.table)
library(stringi)
"%ni%" <- Negate("%in%")

munge_plib <- function(lib_id){
  plib <- fread(paste0("../../../phage_atac_large_data_files/input/phage-evolution-sequencing/data/",lib_id,"_parsed.csv.gz"))
  sdf <- plib[,.N,c("CDR1nt", "CDR2nt", "CDR3nt", "CDR1aa", "CDR2aa", "CDR3aa")]
  sdf2 <- data.frame(
    ID = paste0(sdf[["CDR1nt"]], "_", sdf[["CDR2nt"]], "_", sdf[["CDR3nt"]]),
    aaID = paste0(sdf[["CDR1aa"]], "_", sdf[["CDR2aa"]], "_", sdf[["CDR3aa"]]),
    N = sdf[["N"]]
  ) %>% mutate(perc = N/sum(N)*100) %>%
    arrange(desc(N))
  colnames(sdf2) <- c("ID", "aaID", paste0(lib_id, "_N"), paste0(lib_id, "_perc"))
  sdf2
}
pl3 <- munge_plib("PLIB3") 
pl2 <- munge_plib("PLIB2") 
pl1 <- munge_plib("PLIB1") 
pl0 <- munge_plib("PLIB0") 


process_clone_library <- function(nBC){
  guess <- ceiling(log10(dim(nBC)[1]))
  nBC_keep <- nBC; nBC_keep$cloneID <- ""; nBC_keep$nNew <- 0
  
  #
  library(stringdist)
  h_threshold <- 2
  # Loop through and eat up barcodes
  idx <- 1
  while(idx < 1001){
    barcode <- as.character(nBC[1,1])
    d <- stringdist(barcode, nBC[[1]], method = "hamming")
    boo <- d <= h_threshold
    boo_ids <- nBC$ID[boo]
    clone_barcode <- dropBarcode <- paste0("clone", formatC(idx, width=guess, flag="0", digits = 20), "_N", sprintf("%09d", sum(boo)))
    nBC_keep[nBC_keep$ID %in% boo_ids, "cloneID"] <- clone_barcode
    nBC_keep[nBC_keep$ID %in% boo_ids, "nNew"] <- c(sum(nBC_keep[nBC_keep$ID %in% boo_ids, 3]), 
                                                    rep(0, length(boo_ids) -1))
    idx <- idx + 1
    
    # Remove barcodes that we've dealt with
    nBC <- nBC[nBC$ID %ni% boo_ids,]
  }
  nBC_keep
  
}

collapsed_PLIB3 <- process_clone_library(pl3)
collapsed_PLIB2 <- process_clone_library(pl2)
collapsed_PLIB1 <- process_clone_library(pl1)
collapsed_PLIB0 <- process_clone_library(pl0)

sum(head(collapsed_PLIB3$nNew, 500))/sum(collapsed_PLIB3$PLIB3_N)*100
sum(head(collapsed_PLIB2$nNew, 500))/sum(collapsed_PLIB2$PLIB2_N)*100
sum(head(collapsed_PLIB1$nNew, 500))/sum(collapsed_PLIB1$PLIB1_N)*100
sum(head(collapsed_PLIB0$nNew, 500))/sum(collapsed_PLIB0$PLIB0_N)*100

sum(head(collapsed_PLIB3$nNew, 1000))/sum(collapsed_PLIB3$PLIB3_N)*100
sum(head(collapsed_PLIB2$nNew, 1000))/sum(collapsed_PLIB2$PLIB2_N)*100
sum(head(collapsed_PLIB1$nNew, 1000))/sum(collapsed_PLIB1$PLIB1_N)*100
sum(head(collapsed_PLIB0$nNew, 1000))/sum(collapsed_PLIB0$PLIB0_N)*100


top1000 <- c(0.1459182,9.434512,44.98175,72.63229)
top500 <- c(0.09380088,6.294692,35.25877, 64.38532)
data.frame(
  what = c("zIn", "zIn", "zIn", "zIn", "yHigh", "yHigh", "yHigh","yHigh", "Out","Out",  "Out", "Out"),
  round = c(0:3, 0:3, 0:3),
  value = c(top500, top1000-top500,100-( top1000))
) %>%
  ggplot(aes(x = round, y = value, fill = what)) +
  geom_bar(stat = "identity") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("lightgrey", "dodgerblue3", "firebrick")) -> outp
outp

cowplot::ggsave2(outp, file = "../output/stacked_bar_library-collapsed.pdf", width = 2.7, height = 2)



# Write percentages to xlsx
library(xlsx)
fn <- paste0("../output/topCollapsedClones_percent-top1k.xlsx")
write.xlsx(collapsed_PLIB0 %>% filter(nNew > 0) %>% head(1000) , file=fn, sheetName="PLIB0", row.names=TRUE)
write.xlsx(collapsed_PLIB1 %>% filter(nNew > 0) %>% head(1000) , file=fn, sheetName="PLIB1", append=TRUE, row.names=TRUE)
write.xlsx(collapsed_PLIB2 %>% filter(nNew > 0) %>% head(1000) , file=fn, sheetName="PLIB2", append=TRUE, row.names=TRUE)
write.xlsx(collapsed_PLIB3 %>% filter(nNew > 0) %>% head(1000) , file=fn, sheetName="PLIB3", append=TRUE, row.names=TRUE)

