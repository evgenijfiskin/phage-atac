library(Signac)
library(Seurat)
library(dplyr)
library(BuenColors)

load("../../../phage_atac_large_data_files/output/asap_phage_comparison/15Feb2021_analyzed_seurat_objects.rda")

CD4anno_cluster <- c(0,1,8)
CD8anno_cluster <- c(5,6)

# Compute the PR curves
ddt <- data.frame(
  coembed@meta.data,
  t(coembed@assays$PDT@scale.data),
  t(coembed@assays$ADT@scale.data),
  t(coembed@assays$ASAP@scale.data),
  CD4_cluster = coembed@meta.data$seurat_clusters %in% CD4anno_cluster,
  CD8_cluster = coembed@meta.data$seurat_clusters %in% CD8anno_cluster
)

cite <- ddt %>% dplyr::filter(tech4 == "CITE")
asap <- ddt %>% dplyr::filter(tech4 == "ASAP")
phage <- ddt %>% dplyr::filter(tech4 == "PHAGE")
phage_and_asap <- ddt %>% dplyr::filter(tech4 == "ASAPandPHAGE")

library(precrec)
CITE_CD4 <- evalmod(scores = cite$CD4, labels = cite$CD4_cluster)
CITE_CD8 <- evalmod(scores = cite$CD8, labels = cite$CD8_cluster)
ASAP_CD4 <- evalmod(scores = asap$asapCD4, labels = asap$CD4_cluster)
ASAP_CD8 <- evalmod(scores = asap$asapCD8, labels = asap$CD8_cluster)
PHAGE_CD4 <- evalmod(scores = phage$NB17, labels = phage$CD4_cluster)
PHAGE_CD8 <- evalmod(scores = phage$NB25, labels = phage$CD8_cluster)

coASAP_CD4 <- evalmod(scores = phage_and_asap$asapCD4, labels = phage_and_asap$CD4_cluster)
coASAP_CD8 <- evalmod(scores = phage_and_asap$asapCD8, labels = phage_and_asap$CD8_cluster)
coPHAGE_CD4 <- evalmod(scores = phage_and_asap$NB17, labels = phage_and_asap$CD4_cluster)
coPHAGE_CD8 <- evalmod(scores = phage_and_asap$NB25, labels = phage_and_asap$CD8_cluster)


# ROC plots
cd4_df <- rbind(data.frame((CITE_CD4[[1]][[1]][c(1,2)]),tech = "CITE", marker = "CD4"),
      data.frame((ASAP_CD4[[1]][[1]][c(1,2)]),tech = "ASAP", marker = "CD4"),
      data.frame((PHAGE_CD4[[1]][[1]][c(1,2)]),tech = "Phage", marker = "CD4"),
      data.frame((coPHAGE_CD4[[1]][[1]][c(1,2)]),tech = "coPhage", marker = "CD4")
)

p1 <- ggplot(cd4_df, aes(x = x, y = y, color = tech)) +
  geom_line() + labs(x = "1-Specificity", y = "Sensitivity") +
  scale_color_manual(values = c("Phage"="#ec1c24", "ASAP"="#314f9b", "coPhage"="#8362a5","CITE" = "#000000")) +
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/ROCcurve_CD4.pdf", width = 2, height = 2)

cd8_df <- rbind(data.frame((CITE_CD8[[1]][[1]][c(1,2)]),tech = "CITE", marker = "CD8"),
                data.frame((ASAP_CD8[[1]][[1]][c(1,2)]),tech = "ASAP", marker = "CD8"),
                data.frame((PHAGE_CD8[[1]][[1]][c(1,2)]),tech = "Phage", marker = "CD8"),
                data.frame((coPHAGE_CD8[[1]][[1]][c(1,2)]),tech = "coPhage", marker = "CD8")
)

p2 <- ggplot(cd8_df, aes(x = x, y = y, color = tech)) +
  geom_line() + labs(x = "1-Specificity", y = "Sensitivity") +
  scale_color_manual(values = c("Phage"="#ec1c24", "ASAP"="#314f9b", "coPhage"="#8362a5","CITE" = "#000000")) +
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none")
cowplot::ggsave2(p2, file = "../plots/ROCcurve_CD8.pdf", width = 2, height = 2)


# Summary bar graph
results_df <- rbind(
  data.frame(tech = "CITE", marker = "CD4", auc(CITE_CD4)[,c(3,4)]),
  data.frame(tech = "CITE", marker = "CD8", auc(CITE_CD8)[,c(3,4)]),
  data.frame(tech = "ASAP", marker = "CD4", auc(ASAP_CD4)[,c(3,4)]),
  data.frame(tech = "ASAP", marker = "CD8", auc(ASAP_CD8)[,c(3,4)]),
  data.frame(tech = "PHAGE", marker = "CD4", auc(PHAGE_CD4)[,c(3,4)]),
  data.frame(tech = "PHAGE", marker = "CD8", auc(PHAGE_CD8)[,c(3,4)]),
 # data.frame(tech = "coASAP", marker = "CD4", auc(coASAP_CD4)[,c(3,4)]),
#  data.frame(tech = "coASAP", marker = "CD8", auc(coASAP_CD8)[,c(3,4)]),
  data.frame(tech = "coPHAGE", marker = "CD4", auc(coPHAGE_CD4)[,c(3,4)]),
  data.frame(tech = "coPHAGE", marker = "CD8", auc(coPHAGE_CD8)[,c(3,4)])
)
results_df %>% dplyr::filter(curvetypes == "ROC") %>%
  mutate(aucor = round(aucs, 2)) %>% 
  ggplot(aes(x = marker, y = aucs, fill = tech)) +
  geom_bar(stat = "identity", position="dodge", color = "black") +
  labs(x = "", y = "AUROC") +
  geom_text(aes(label = aucor), vjust = -0.5,  position = position_dodge(0.9), size = 2) +
  scale_y_continuous() + pretty_plot() + L_border() + theme(legend.position = "bottom") +
  scale_fill_manual(values = jdb_palette("corona")) -> pOut
cowplot::ggsave2(pOut, file = "../plots/ROC_barchart.pdf", width = 3, height = 3.5)
