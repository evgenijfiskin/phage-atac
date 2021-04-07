# Load libraries
library(dplyr)
library(data.table)

dd <- fread("../data/commercial/Fig1_singlecell.csv.gz") %>% filter(cell_id != "None")
mean(dd$DNase_sensitive_region_fragments/dd$passed_filters)
summary(dd$passed_filters)
