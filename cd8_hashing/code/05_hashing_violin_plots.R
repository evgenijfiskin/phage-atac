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
CD8counts <- read.csv(file = "../data/EF_CD8_hashing_analysis.csv")

CD8counts$PDT_classification.global <- factor(CD8counts$PDT_classification.global, levels= c("Doublet","Negative","Singlet"))

Fig2e <- ggplot(CD8counts, aes(x=PDT_classification.global, y=passed_filters, fill=PDT_classification.global)) + ylim(0,20000) +
  geom_violin(trim=FALSE, show.legend = FALSE)+ theme_classic() + scale_fill_manual(values=c("#EE322A", "#D3D3D3", "#4E65AF")) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median,
                                                                                                                                              geom = "crossbar", width = 0.7, color="black",show.legend = FALSE) + geom_jitter(shape=16, size=0.1, position=position_jitter(0.2),show.legend = FALSE)
Fig2f <- ggplot(CD8counts, aes(x=PDT_classification.global, y=nCount_PDT, fill=PDT_classification.global)) + ylim(0,20000) +
  geom_violin(trim=FALSE, show.legend = FALSE)+ theme_classic() + scale_fill_manual(values=c("#EE322A", "#D3D3D3", "#4E65AF")) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median,
                                                                                                                                              geom = "crossbar", width = 0.7, color="black",show.legend = FALSE) + geom_jitter(shape=16, size=0.1, position=position_jitter(0.2),show.legend = FALSE)
EDFig11e <- ggplot(CD8counts, aes(x=CD8pos, y=maxCD8, fill=CD8pos)) +
  geom_violin(trim=FALSE, show.legend = FALSE)+ theme_classic() + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median,
                                                                                                                                              geom = "crossbar", width = 0.7, color="black",show.legend = FALSE) + geom_jitter(shape=16, size=0.1, position=position_jitter(0.2),show.legend = FALSE)
EDFig11f <- ggplot(CD8counts, aes(x=CD8pos, y=log10(passed_filters), fill=CD8pos)) + ylim(0,6) +
  geom_violin(trim=FALSE, show.legend = FALSE)+ theme_classic() + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median,
                                                                               geom = "crossbar", width = 0.7, color="black",show.legend = FALSE) + geom_jitter(shape=16, size=0.1, position=position_jitter(0.2),show.legend = FALSE)