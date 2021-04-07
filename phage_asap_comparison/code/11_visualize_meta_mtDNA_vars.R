library(BuenColors)
source("../../global_functions/variant_enrichment.R")

ref_all <- fread("../../global_functions/references/chrM_refAllele.txt")
vars_clone <- readRDS("../output/Phage_ASAP_all_possibleVars_meta.rds")
prop_df <- get_enrich_mutation_df(rownames(vars_clone), ref_all)

# Visualize the nucleotide bias
p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot(fontsize = 8) + L_border() + 
  theme(axis.title.x=element_blank(),
        axis.text.x =element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide (n = 518 mutations)", y = "Substitution Rate (Expected / Observed)")
cowplot::ggsave2(p1, file = "../plots/Phage_n518_mito_signature.pdf", width = 4, height = 2.4)




misc_df <- 
  low_freq_df <- misc_df %>% filter(n_cells_conf_detected > 0 & n_cells_conf_detected <= 2 & n_cells_over_20 > 0) %>% filter(n_cells_over_5 <= 5) 
dim(low_freq_df)

p1 <- ggplot(misc_df %>%  filter(n_cells_conf_detected >= 2 ), aes(x = strand_correlation, y = log10(vmr), color = log10(vmr) > -2 & strand_correlation > 0.55)) +
  geom_point(size = 0.4) + scale_color_manual(values = c("black", "firebrick")) +
  labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "bottom") + 
  geom_vline(xintercept = 0.55, linetype = 2) +
  geom_hline(yintercept = -2, linetype = 2) + theme(legend.position = "none")
p1
cowplot::ggsave2(p1, file = "../plots/PhageATAC_var_call_mtscatac.pdf", width = 2, height = 2)