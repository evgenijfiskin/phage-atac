library(BuenColors)
library(data.table)
library(dplyr)

facs <- fread("../data/Intracellular_FACS_data.txt")
facs$classification <- ifelse(facs[[1]] > 5e4, "high", "low")
# Make Plots
pFacs <- ggplot(facs, aes(x = log10(FITC), fill = classification)) +
  geom_histogram(data = facs %>% filter(classification == "low"), aes(y = ..density..), alpha = 0.3, bins = 20) +
  geom_histogram(data = facs %>% filter(classification == "high"), aes(y = ..density..), alpha = 0.3, bins = 20) +
  geom_density(aes(color = classification), fill = NA) +
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(3,7)) +
  scale_color_manual(values = c("#4eabd9","#4e66ae")) +
  scale_fill_manual(values = c("#4eabd9","#4e66ae")) +
  ggtitle("FACS") +  labs(x = "log10 eGFP (FACS)") +
  theme(legend.position = "none") + theme(plot.title = element_text(size = 8))
pFacs

import_phage_counts <- function(){
  dt <- fread(paste0("../data/PA_GFP_21Feb2021.kallisto.counts.tsv"))
  dm <- data.matrix(dt[,c(-1,-3)])
  rownames(dm) <- paste0(dt[[1]],"-1")
  return(t(dm))
}

singlecell <- fread('../data/singlecell.csv', header = TRUE) %>% dplyr::filter(cell_id != "None")
pdt <- import_phage_counts()[6,]
pdt <- pdt[names(pdt) %in% singlecell$barcode]
phage  <- data.frame(
  barcode = names(pdt), 
  intra_eGFP_counts = unname(pdt)
) %>%  mutate(classification = ifelse(intra_eGFP_counts > 100, "high", "low"))

qplot(phage[[2]], bins = 30) + scale_x_log10() + geom_vline(xintercept = 100)
pPhage <- ggplot(phage , aes(x = log10(intra_eGFP_counts), fill = classification)) +
  geom_histogram(data = phage %>% filter(classification == "low"), aes(y = ..density..), alpha = 0.3, bins = 20) +
  geom_histogram(data = phage %>% filter(classification == "high"), aes(y = ..density..), alpha = 0.3, bins = 20) +
  geom_density(aes(color = classification), fill = NA, adjust = 2) +
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("#4eabd9","#4e66ae")) +
  scale_fill_manual(values =  c("#4eabd9","#4e66ae")) +
  ggtitle("Phage") + labs(x = "log10 eGFP (reads)") +
  theme(legend.position = "none") + theme(plot.title = element_text(size = 8))
pPhage
table(facs$classification)
table(phage$classification)

cowplot::ggsave2(cowplot::plot_grid(pFacs, pPhage, ncol = 2), 
                 file = "../plots/histograms_intragfp2-sbs.pdf", width = 3.7, height = 2)
