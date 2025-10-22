library(Seurat)
library(Signac)
library(SeuratWrappers)
library(RENIN)
pbmc <- readRDS("scENGRN/method/pbmc/data/pbmcSCTRENIN_motif.rds")

file_path <- "scENGRN/method/pbmc/data/highly_variable_genes.txt"
highly_variable_genes <- readLines(file_path) 

expr_mat <- GetAssayData(pbmc, assay = "SCT", layer = "data") 
expr_mat <- t(expr_mat)
peak_WNN <-readRDS("scENGRN/method/pbmc/data/weighted_atac.rds")
peak_WNN <- t(peak_WNN)
peak_results <- run_peak_gene(pbmc, expr_mat, peak_WNN, highly_variable_genes, lambda2 = 0.25, max_distance = 5e+05)


aen_lists <- make_aen_lists(peak_results)
gene_names <- names(aen_lists)
spatial_lag_matrix <- readRDS("scENGRN/method/pbmc/data/weighted_rna.rds")
spatial_lag_matrix <- spatial_lag_matrix[gene_names, ]
spatial_lag_matrix <- t(spatial_lag_matrix)
tf_results <- run_tf_gene(pbmc, expr_mat, peak_results, spatial_lag_matrix, gene_names, lambda2 = 0.5,layer = "data")