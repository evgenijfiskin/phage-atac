library(data.table)
library(dplyr)
library(BuenColors)
library(circlize)
library(ComplexHeatmap)
source('variant_calling.R')

# Import output file from mgatk
m <- readRDS("../../../phage_atac_large_data_files/input/CD8_hashing/phage_atac_cd8_mgatk/final/phage_atac_cd8_mgatk.rds")
SE <- m[,m$depth > 5]
mut_se <- call_mutations_mgatk(SE)
misc_df <- data.frame(rowData(mut_se))
filter_df <- misc_df %>%  filter(n_cells_conf_detected >= 10 &  n_cells_over_20 >= 10 & strand_correlation > 0.6 )
dim(filter_df)

# Compare with Gokcen calls
gc <- fread("../data/Gokcen_clusters.tsv.gz")
gc$nu <- ifelse(gc$assignment %in% c("0", "1", "2", "3"), paste0("d", gc$assignment), "doublet")
vec <- gc$nu; names(vec) <- gc$barcode

color_vec <- c("dodgerblue2", "purple3", "red", "green2", "black"); names(color_vec) <- c("d0", "d1", "d2", "d3", "doublet")
ha_col <- HeatmapAnnotation(df = data.frame(Gokcen = vec[colnames(mut_se)]),
                            col = list(Gokcen = color_vec))

# Visualize heatmap
pdf(paste0("../plots/mito_heatmap_GC_clusters_split.pdf"), width=10, height=10)
hm <- Heatmap(data.matrix(assays(mut_se)[["allele_frequency"]][filter_df$variant,]), 
              column_split = vec[colnames(mut_se)],
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              cluster_columns = TRUE,
              top_annotation=ha_col,
              name = "AF",
              cluster_rows = TRUE, 
              show_column_names = FALSE)
hm
dev.off()
