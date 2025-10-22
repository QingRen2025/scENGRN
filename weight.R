# 加载必要的库
library(SummarizedExperiment)
library(Matrix)
library(GenomicRanges)

wnn_matrix <- readRDS("scENGRN/method/pbmc/data/wnn_matrix.rds")

# 对每一行进行标准化（使每一行的和为 1）
wnn_matrix_normalized <- wnn_matrix / rowSums(wnn_matrix)
# 保存标准化后的 WNN 矩阵为 RDS 文件
saveRDS(wnn_matrix_normalized, file = "scENGRN/method/pbmc/data/wnn_matrix_normalized.rds")
# 加载 RDS 文件
wnn_matrix_normalized <- readRDS("scENGRN/method/pbmc/data/wnn_matrix_normalized.rds")

pbmc <- readRDS("scENGRN/method/pbmc/data/pbmc_seurat_object.rds")

rna <- pbmc@assays$SCT@data
# RNA 加权平均
weighted_rna <- rna %*% wnn_matrix_normalized
# 查看结果的维度
dim(weighted_rna)
# 保存 weighted_rna 为 RDS 文件
saveRDS(weighted_rna, file = "scENGRN/method/pbmc/data/weighted_rna.rds")

atac <- pbmc@assays[["ATAC"]]@data
# ATAC 加权平均
weighted_atac <- atac %*% wnn_matrix_normalized
# 查看结果的维度
dim(weighted_atac)
saveRDS(weighted_atac, file = "scENGRN/method/pbmc/data/weighted_atac.rds")


