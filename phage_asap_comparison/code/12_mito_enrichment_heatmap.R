library(data.table)
library(dplyr)
library(BuenColors)
process_library <- function(lib){
  dat <- fread(paste0("../data/Lib",lib,"/EF_PBMC_",lib,"_ATAC_hg38_v12-mtMask.singlecell.csv.gz"))
  dat$library <- lib
  dat %>% filter(passed_filters > 500 & barcode != "NO_BARCODE") -> dat2
  dat2$pct_mito <- dat2$mitochondrial/dat2$total *100
  dat2
}
dt <- rbind(process_library("A"), process_library("B"), process_library("C"))
dt$density <- get_density(log10(dt$passed_filters), dt$pct_mito)

ggplot(dt %>% arrange((density)), aes(x = log10(passed_filters), y = pct_mito, color = density)) +
  geom_point() + scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  facet_wrap(~library) + pretty_plot() + L_border() +
  theme(legend.position = "none") + labs(y = "% Mito", x = "log10 chromatin fragments")
