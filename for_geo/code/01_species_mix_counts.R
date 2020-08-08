library(data.table)

#---
# Commerical 
#--- 

tag_df <- rbind(
  fread("../../species_mix/data/Fig1_phcellbarcodeCOUNTSfiltered_human_clean.csv.gz"),
  fread("../../species_mix/data/Fig1_phcellbarcodeCOUNTSfiltered_mouse_clean.csv.gz"),
  fread("../../species_mix/data/Fig1_phcellbarcodeCOUNTSfiltered_noncells_clean.csv.gz")
)
tag_df$barcode <- paste0(tag_df$V2, "-1")
o_df <- tag_df[,c(3,1)]; colnames(o_df) <- c("barcode", "eGFP_counts")

# Merge with single cell data
tenx_df <- fread("../../species_mix/data/Fig1_singlecell.csv.gz")
merged_df <- left_join(tenx_df, o_df, by = "barcode")
merged_df$eGFP_counts <- ifelse(!is.na(merged_df$eGFP_counts), merged_df$eGFP_counts, 0)
write.table(merged_df, file = "../out/SpeciesMix_eGFP_Commercial_singlecell_counts.csv", 
            sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

#---
# Evovled 
#--- 

tag_df <- rbind(
  fread("../../species_mix/data/Fig2_C5_barcodeCOUNTSfiltered_human_clean.csv.gz"),
  fread("../../species_mix/data/Fig2_C5_barcodeCOUNTSfiltered_mouse_clean.csv.gz"),
  fread("../../species_mix/data/Fig2_C5_barcodeCOUNTSfiltered_noncells_clean.csv.gz")
)
tag_df$barcode <- paste0(tag_df$V2, "-1")
o_df <- tag_df[,c(3,1)]; colnames(o_df) <- c("barcode", "eGFP_counts")

# Merge with single cell data
tenx_df <- fread("../../species_mix/data/Fig2_singlecell.csv.gz")
merged_df <- left_join(tenx_df, o_df, by = "barcode")
merged_df$eGFP_counts <- ifelse(!is.na(merged_df$eGFP_counts), merged_df$eGFP_counts, 0)
write.table(merged_df, file = "../out/SpeciesMix_eGFP_Evolved_singlecell_counts.csv", 
            sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
