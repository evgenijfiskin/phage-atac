library(BuenColors)
library(dplyr)
library(data.table)


munge_plib <- function(lib_id){
  plib <- fread(paste0("../../../phage_atac_large_data_files/input/phage-evolution-sequencing/data/",lib_id,"_parsed.csv.gz"))
  sdf <- plib[,.N,c("CDR1aa", "CDR2aa", "CDR3aa")]
  sdf2 <- data.frame(
    ID = paste0(sdf[["CDR1aa"]], "_", sdf[["CDR2aa"]], "_", sdf[["CDR3aa"]]),
    N = sdf[["N"]]
  ) %>% mutate(perc = N/sum(N)*100) %>%
    arrange(desc(N))
  colnames(sdf2) <- c("ID", paste0(lib_id, "_N"), paste0(lib_id, "_perc"))
  sdf2
}
pl3 <- munge_plib("PLIB3") %>% mutate(clone_ID = paste0("clone", as.character(1:n())))
pl2 <- munge_plib("PLIB2") 
pl1 <- munge_plib("PLIB1") 
pl0 <- munge_plib("PLIB0") 

m_all_df <- merge(pl3, pl2, by = "ID", all.x = TRUE) %>%
  merge(pl1, by = "ID", all.x = TRUE) %>%
  merge(pl0, by = "ID", all.x = TRUE) 
all_df <- m_all_df %>% arrange(desc(PLIB3_perc))

write.table(head(all_df,500), file = "../output/top500clones_quant.csv", 
            sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

abundance_rank_df <- data.frame(
  rank = c(1:dim(pl0)[1], 1:dim(pl1)[1], 1:dim(pl2)[1], 1:dim(pl3)[1]),
  abundance = c(pl0[[2]], pl1[[2]], pl2[[2]], pl3[[2]]),
  what = c(rep("p0", dim(pl0)[1]), rep("p1", dim(pl1)[1]), rep("p2", dim(pl2)[1]), rep("p3", dim(pl3)[1]))
)


ggplot(abundance_rank_df %>% dplyr::filter(rank < 100000), aes(x = rank, y = abundance, color = what)) + 
  geom_point() + scale_x_log10() + scale_y_log10()

# Stacked bar
ss_df <- head(all_df,500)
sum(ss_df$PLIB3_perc, na.rm = TRUE)
sum(ss_df$PLIB2_perc, na.rm = TRUE)
sum(ss_df$PLIB1_perc, na.rm = TRUE)
sum(ss_df$PLIB0_perc, na.rm = TRUE)
vecvalue <- c( 0.02789233, 3.635823, 25.73607, 54.88258)

m_all_df %>% arrange(desc(PLIB0_perc)) %>% head(500) %>% pull(PLIB0_perc) %>% sum()
m_all_df %>% arrange(desc(PLIB1_perc)) %>% head(500) %>% pull(PLIB1_perc) %>% sum()
m_all_df %>% arrange(desc(PLIB2_perc)) %>% head(500) %>% pull(PLIB2_perc) %>% sum()
m_all_df %>% arrange(desc(PLIB3_perc)) %>% head(500) %>% pull(PLIB3_perc) %>% sum()

data.frame(
  what = c("zIn", "zIn", "zIn", "zIn", "Out", "Out", "Out", "Out"),
  round = c(0:3, 0:3),
  value = c(vecvalue, 100-vecvalue)
) %>%
  ggplot(aes(x = round, y = value, fill = what)) +
  geom_bar(stat = "identity", color = "black") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("lightgrey", "firebrick")) -> outp
cowplot::ggsave2(outp, file = "../output/stacked_bar_library.pdf", width = 2.7, height = 2)

