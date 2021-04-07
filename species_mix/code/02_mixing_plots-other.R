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
EGFPcounts <- read.csv(file = "../data/EF_293T_3T3_mixing_analysis.csv")

EGFPcounts$Type <- factor(EGFPcounts$Type, levels= c("Human EGFP+","Human EGFP-","Mouse"))

Fig1g <- ggplot(EGFPcounts, aes(x=Type, y=EGFPlog, fill=Type)) + ylim(0,5.5) +
  geom_violin(trim=FALSE, show.legend = FALSE)+ theme_classic() + scale_fill_manual(values=c("#429DD7", "#4E65AF", "#EE322A")) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median,
                                                                                                                                              geom = "crossbar", width = 0.7, color="black",show.legend = FALSE) + geom_jitter(shape=16, size=0.1, position=position_jitter(0.2),show.legend = FALSE)

Fig1h <- ggplot(EGFPcounts, aes(x=Type, y=Totallog, fill=Type)) + ylim(0,5.5) +
  geom_violin(trim=FALSE, show.legend = FALSE)+ theme_classic() + scale_fill_manual(values=c("#429DD7", "#4E65AF", "#EE322A")) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median,
                                                                                                                                              geom = "crossbar", width = 0.7, color="black",show.legend = FALSE) + geom_jitter(shape=16, size=0.1, position=position_jitter(0.2),show.legend = FALSE)

EDFig7b <- ggplot(EGFPcounts, aes(x=Type, y=Percentage.in.peaks, fill=Type)) + ylim(0,1) +
  geom_violin(trim=FALSE, show.legend = FALSE)+ theme_classic() + scale_fill_manual(values=c("#429DD7", "#4E65AF", "#EE322A")) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median,
                                                                                                                                              geom = "crossbar", width = 0.7, color="black",show.legend = FALSE) + geom_jitter(shape=16, size=0.1, position=position_jitter(0.2),show.legend = FALSE)

EDFig7c <- ggplot(EGFPcounts, aes(x=Type, y=Percentage.overlapping.TSS, fill=Type)) + ylim(0,1) +
  geom_violin(trim=FALSE, show.legend = FALSE)+ theme_classic() + scale_fill_manual(values=c("#429DD7", "#4E65AF", "#EE322A")) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median,
                                                                                                                                              geom = "crossbar", width = 0.7, color="black",show.legend = FALSE) + geom_jitter(shape=16, size=0.1, position=position_jitter(0.2),show.legend = FALSE)