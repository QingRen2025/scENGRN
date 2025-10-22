library(SummarizedExperiment)
library(Matrix)
library(GenomicRanges)

wnn_matrix <- readRDS("wnn_matrix.rds")
# Normalize each row (so that the sum of each row equals 1)
wnn_matrix_normalized <- wnn_matrix / rowSums(wnn_matrix)

pbmc <- readRDS("pbmc_seurat_object.rds")

rna <- pbmc@assays$SCT@data
# RNA weighted average
weighted_rna <- rna %*% wnn_matrix_normalized
saveRDS(weighted_rna, file = "weighted_rna.rds")

atac <- pbmc@assays[["ATAC"]]@data
# ATAC weighted average
weighted_atac <- atac %*% wnn_matrix_normalized
saveRDS(weighted_atac, file = "weighted_atac.rds")


