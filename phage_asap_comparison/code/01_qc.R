library(dplyr)
library(data.table)
library(BuenColors)
library(ggrastr)
library(SummarizedExperiment)
library(Matrix)

source("../../global_functions/estimateLibraryComplexity.R")

set.seed(1)
# Get top n cell qc metrics per run
importQC <- function(raw_name, name, n = 1000){
  
  dt <- fread(raw_name) %>% 
    dplyr::filter(cell_id != "None") %>% 
    mutate(Experiment = name) %>%
    mutate(DNaseProp = DNase_sensitive_region_fragments/passed_filters) %>%
    mutate(MitoProp = mitochondrial/total) %>% 
    mutate(duplicateProp = duplicate / total) %>%
    mutate(TSSProp = TSS_fragments/passed_filters) %>% data.frame()
  dt$LibraryComplexity <- sapply(1:dim(dt)[1], function(i){
    estimateLibrarySize(dt[i,"total"],dt[i,"passed_filters"])})
  
  top_nk <- dt[sample(1:dim(dt)[1], n),]
  top_nk$barcode <- as.character(top_nk$barcode)
  return(top_nk)
  
}

# Import the original screen conditions into one data frame
screen_df <- rbindlist(list(
  importQC("../data/LibA/EF_PBMC_A_ATAC_hg38_v12-mtMask.singlecell.csv.gz", "ASAP"),
  importQC("../data/LibB/EF_PBMC_B_ATAC_hg38_v12-mtMask.singlecell.csv.gz", "CoStain"),
  importQC("../data/LibC/EF_PBMC_C_ATAC_hg38_v12-mtMask.singlecell.csv.gz", "Phage")
)) %>% data.frame()

mwtest <- function(exp1, exp2, attribute){
  v1 <- screen_df[screen_df$Experiment == exp1,attribute]
  v2 <- screen_df[screen_df$Experiment == exp2,attribute]
  
  data.frame(exp1, exp2, attribute, pvalue = (wilcox.test(v1, v2)$p.value))
}
mwtest("ASAP", "CoStain", "LibraryComplexity")
mwtest("ASAP", "Phage", "LibraryComplexity")
mwtest("Phage", "CoStain", "LibraryComplexity")

mwtest("ASAP", "CoStain", "DNaseProp")
mwtest("ASAP", "Phage", "DNaseProp")
mwtest("Phage", "CoStain", "DNaseProp")

mwtest("ASAP", "CoStain", "TSSProp")
mwtest("ASAP", "Phage", "TSSProp")
mwtest("Phage", "CoStain", "TSSProp")

mwtest("ASAP", "CoStain", "MitoProp")
mwtest("ASAP", "Phage", "MitoProp")
mwtest("Phage", "CoStain", "MitoProp")


# Summarize values for the initial screen of conditions
screen_df %>% group_by(Experiment) %>%
  summarize(PCT_Mito = median(MitoProp)*100, 
            PCT_DNase = median(DNaseProp)*100, 
            mean(DNaseProp)*100,
            PCT_TSS = median(TSSProp)*100,
            Chromatin_Complexity = median(LibraryComplexity), count =n ())

# Visualize essential attributes -- mito and accessible chromatin enrichment
pA <- ggplot(screen_df, aes(x = Experiment, y = MitoProp*100)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.2) +
  coord_cartesian(ylim = c(0, 70)) + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 7) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "% mtDNA fragments") +
  scale_color_manual(values = "black") + L_border() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

pB <- ggplot(screen_df, aes(x = Experiment, y = DNaseProp*100)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.2) +
  coord_cartesian(ylim = c(25, 95)) + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 7) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "% reads in DNase") +
  scale_color_manual(values = "black") + L_border() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

cowplot::ggsave2(cowplot::plot_grid(pA, pB, nrow = 1), file = "../plots/Phage-ASAP_Panel1BC.pdf", width = 3.7, height = 2)

# Summarize additional plots for the supplement
pC <- ggplot(screen_df, aes(x = Experiment, y = TSSProp*100)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.2) + 
  coord_cartesian(ylim = c(0, 75)) + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 7) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "% reads overlap TSS") +
  scale_color_manual(values = "black") + L_border() + theme(axis.text.x = element_text(angle = 0, hjust = 1))

cowplot::ggsave2(pC, file = "../plots/Phage-ASAP_tss_rates.pdf", width = 1.8, height = 1.8)

pD <-  ggplot(screen_df, aes(x = Experiment, y = log10(LibraryComplexity))) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.2) + 
  coord_cartesian(ylim = c(3, 5.5))  + scale_y_continuous(expand = c(0,0)) + 
  pretty_plot(fontsize = 6) + theme(legend.position = "bottom") + labs(x = "", color = "", y = "log10 library size") +
  scale_color_manual(values = "black") + L_border()+ theme(axis.text.x = element_text(angle = 0, hjust = 1))

cowplot::ggsave2(pD, file = "../plots/Phage-ASAP_qc_library_complexity.pdf", width = 1.8, height = 1.8)



