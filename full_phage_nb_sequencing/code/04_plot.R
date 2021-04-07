library(BuenColors)
library(dplyr)

t500 <- read.csv("../output/top500clones_quant.csv")

reshape2::melt(t500[,c("clone_ID","PLIB0_perc", "PLIB1_perc", "PLIB2_perc", "PLIB3_perc")], id.vars = "clone_ID") %>%
  filter(clone_ID %in% paste0("clone", 1:500)) -> raw_df

plot_df <- data.frame(
  lib = gsub("_perc", "", raw_df$variable ),
  clone_ID =  raw_df$clone_ID,
  value = ifelse(is.na(raw_df$value), 0.000001, raw_df$value)
)
  
ggplot(plot_df, aes(x = lib, y = value, color = clone_ID, group = clone_ID)) +
  geom_point() + geom_line() +
  pretty_plot() + L_border() +
  labs(x = "Library", y = "% abundance in library") +
  theme(legend.position = "none") 
