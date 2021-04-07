library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(data.table)
library(BuenColors)
"%ni%" <- Negate("%in%")
call_mutations_mgatk_expanded <- function(SE, stabilize_variance = TRUE, low_coverage_threshold = 10){
  # Determinie key coverage statistics every which way
  cov <- assays(SE)[["coverage"]]
  ref_allele <- toupper(as.character(rowRanges(SE)$refAllele))
  # Process mutation for one alternate letter
  process_letter <- function(letter){
    print(letter)
    boo <- ref_allele != letter & ref_allele != "N"
    pos <- start(rowRanges(SE))
    variant_name <- paste0(as.character(pos), ref_allele, ">", letter)[boo]
    nucleotide <- paste0(ref_allele, ">", letter)[boo]
    position_filt <- pos[boo]
    # Single cell functions
    getMutMatrix <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    getMutMatrix_fw  <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_fw")]]) / cov_fw)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    getMutMatrix_rev  <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_rev")]]) / cov_rev)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    # Bulk functions
    getBulk <- function(letter){
      vec <- (Matrix::rowSums(assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / Matrix::rowSums(cov))[boo]
      return(vec)
    }
    rowVars <- function(x, ...) {
      Matrix::rowSums((x - Matrix::rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
    }
    update_missing_w_zero <- function(vec){
      ifelse(is.na(vec)  | is.nan(vec), 0, vec)
    }
    # Set up correlation per non-zero mutation based on the strands
    dt <- merge(data.table(Matrix::summary(assays(SE)[[paste0(letter, "_counts_fw")]][boo,])), 
                data.table(Matrix::summary(assays(SE)[[paste0(letter, "_counts_rev")]][boo,])), 
                by.x = c("i", "j"), by.y = c("i", "j"), 
                all = TRUE)[x.x >0 | x.y >0]
    dt$x.x <- update_missing_w_zero(dt$x.x)
    dt$x.y <- update_missing_w_zero(dt$x.y)
    dt2 <- data.table(variant = variant_name[dt[[1]]],
                      cell_idx = dt[[2]], 
                      forward = dt[[3]],
                      reverse = dt[[4]])
    rm(dt)
    cor_dt <- dt2[, .(cor = cor(c(forward), c(reverse), method = "pearson", use = "pairwise.complete")), by = list(variant)]
    # Put in vector for convenience
    cor_vec_val <- cor_dt$cor
    names(cor_vec_val) <- as.character( cor_dt$variant )
    # Compute the single-cell data
    mat <- getMutMatrix(letter)
    mmat <- sparseMatrix(
      i = c(summary(mat)$i,dim(mat)[1]),
      j = c(summary(mat)$j,dim(mat)[2]),
      x = c(update_missing_w_zero(summary(mat)$x), 0)
    )
    rownames(mmat) <- variant_name
    # Get max heteroplasmy per variant
    max_dt <- data.table(summary(mmat))[,.(max = max(x)),by = i]
    max_vec <- max_dt$max; names(max_vec) <- variant_name[as.numeric(max_dt$i)]
    max_vec_all_var <- update_missing_w_zero(max_vec[as.character(rownames(mmat))])
    # Compute bulk statistics
    mean = update_missing_w_zero(getBulk(letter))
    # Stablize variances by replacing low coverage cells with mean
    if(stabilize_variance){
      # Get indices of cell/variants where the coverage is low and pull the mean for that variant
      idx_mat <- which(data.matrix(cov[boo,] < low_coverage_threshold), arr.ind = TRUE)
      idx_mat_mean <- mean[idx_mat[,1]]
      # Now, make sparse matrices for quick conversion
      ones <- 1 - sparseMatrix(
        i = c(idx_mat[,1], dim(mmat)[1]),
        j = c(idx_mat[,2], dim(mmat)[2]),
        x = 1
      )
      means_mat <- sparseMatrix(
        i = c(idx_mat[,1], dim(mmat)[1]),
        j = c(idx_mat[,2], dim(mmat)[2]),
        x = c(idx_mat_mean, 0)
      )
      mmat2 <- mmat * ones + means_mat
      variance = rowVars(mmat2)
      rm(mmat2); rm(ones); rm(means_mat); rm(idx_mat); rm(idx_mat_mean)
    } else {
      variance = rowVars(mmat)
    }
    detected <- (assays(SE)[[paste0(letter, "_counts_fw")]][boo,] >= 2) + (assays(SE)[[paste0(letter, "_counts_rev")]][boo,] >=2 )
    # Compute per-mutation summary statistics
    var_summary_df <- data.frame(
      position = position_filt,
      nucleotide = nucleotide, 
      variant = variant_name,
      vmr = round(variance/(mean + 0.00000000001),8),
      mean = round(mean,7),
      variance = round(variance,7),
      n_cells_conf_detected = Matrix::rowSums(detected == 2),
      n_cells_over_5 = Matrix::rowSums(mmat >= 0.05), 
      n_cells_over_10 = Matrix::rowSums(mmat >= 0.10),
      n_cells_over_20 = Matrix::rowSums(mmat >= 0.20),
      n_cells_over_95 = Matrix::rowSums(mmat >= 0.95),
      max_heteroplasmy = round(max_vec_all_var,3),
      strand_correlation = round(cor_vec_val[variant_name],3),
      mean_coverage = round(Matrix::rowMeans(cov)[boo], 3),
      stringsAsFactors = FALSE, row.names = variant_name
    )
    se_new <- SummarizedExperiment(
      rowData = var_summary_df, 
      colData = colData(SE), 
      assays = list(allele_frequency = mmat, coverage = cov[boo,])
    )
    return(se_new)
  }
  return(SummarizedExperiment::rbind(process_letter("A"), process_letter("C"), process_letter("G"), process_letter("T")))
}
libA <- readRDS("../../../phage_atac_large_data_files/input/asap_comparison/EF_PBMC_A_ATAC_hg38_v12-mtMask_mgatk.rds") 
libB <- readRDS("../../../phage_atac_large_data_files/input/asap_comparison/EF_PBMC_B_ATAC_hg38_v12-mtMask_mgatk.rds")
libC <- readRDS("../../../phage_atac_large_data_files/input/asap_comparison/EF_PBMC_C_ATAC_hg38_v12-mtMask_mgatk.rds")
colnames(libB) <- gsub("-1", "-2", colnames(libB))
colnames(libC) <- gsub("-1", "-3", colnames(libC))

table(libA$depth > 10)
table(libB$depth > 10)
table(libC$depth > 10)

all_mgatk <- call_mutations_mgatk(cbind(
  libA[,libA$depth > 10], libB[,libB$depth > 10], libC[,libC$depth > 10]
))
data.frame(rowData(all_mgatk)) %>%  dplyr::filter(n_cells_conf_detected >= 3 & strand_correlation > 0.65 & log10(vmr) > -2 ) %>% dim()

data.frame(rowData(all_mgatk)) %>%  dplyr::filter(n_cells_conf_detected >= 1 & strand_correlation > 0.65 & log10(vmr) > -2 )-> all_possible

data.frame(rowData(all_mgatk)) %>%  dplyr::filter(n_cells_conf_detected >= 3 & strand_correlation > 0.65 & log10(vmr) > -2 ) %>%
  pull(variant) -> for_plotting

saveRDS(assays(all_mgatk)[["allele_frequency"]][for_plotting,], file = "../output/Phage_ASAP_visualize_coverage.rds")
saveRDS(all_possible, file = "../output/Phage_ASAP_all_possibleVars_meta.rds")
saveRDS(data.frame(rowData(all_mgatk)), file = "../output/Phage_ASAP_everything_meta.rds")
