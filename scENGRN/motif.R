library(Seurat)
library(Signac)
library(SeuratWrappers)
pbmc <- readRDS("pbmc_seurat_object.rds")

library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(chromVARmotifs)

pwm_list <- chromVARmotifs::human_pwms_v2
pbmc <- AddMotifs(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm_list)
saveRDS(pbmc, file = "pbmc_motif.rds")
