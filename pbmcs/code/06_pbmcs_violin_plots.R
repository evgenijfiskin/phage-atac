# Load libraries
library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)
library(scales)
library(cowplot)
library(gplots)
library(GSA)
library(stringr)
library(DropletUtils)
library(devtools)  
library(plyr)

# read csv file
PBMCcounts <- read.csv(file = "../data/EF_PBMC_analysis.csv")

PBMCcounts$Antigen <- factor(PBMCcounts$Antigen, levels= c("CD16","CD4","CD8","EGFP"))

EDFig10e <- ggplot(PBMCcounts, aes(x=Antigen, y=log10(Counts), fill=Antigen)) +
  geom_violin(trim=FALSE, show.legend = FALSE)+ theme_classic() + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median,
                                                                               geom = "crossbar", width = 0.7, color="black",show.legend = FALSE) + geom_jitter(shape=16, size=0.1, position=position_jitter(0.2),show.legend = FALSE)
