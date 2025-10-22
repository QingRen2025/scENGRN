library(Seurat)
library(Signac)
library(SeuratWrappers)
library(RENIN)
pbmc <- readRDS("pbmc_motif.rds")

file_path <- "/highly_variable_genes.txt"
highly_variable_genes <- readLines(file_path) 

RNA <- GetAssayData(pbmc, assay = "SCT", layer = "data") 
RNA <- t(RNA)
peak_WNN <-readRDS("weighted_atac.rds")
peak_WNN <- t(peak_WNN)
peak_gene <- run_peak_gene(pbmc, RNA, peak_WNN, highly_variable_genes, lambda2 = 0.25, max_distance = 5e+05)


aen_lists <- make_aen_lists(peak_gene)
gene_names <- names(aen_lists)
RNA_WNN <- readRDS("weighted_rna.rds")
RNA_WNN <- RNA_WNN[gene_names, ]
RNA_WNN <- t(RNA_WNN)
tf_gene <- run_tf_gene(pbmc, RNA, peak_gene, RNA_WNN, gene_names, lambda2 = 0.5,layer = "data")