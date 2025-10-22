# scENGRN
This pipeline utilizes elastic net regression to infer gene regulatory networks from single-cell RNA-seq and ATAC-seq data. The pipeline is divided into the following steps:
## 1. Raw Data Preparation
Paired gene expression data and chromatin accessibility data are used as input.

## 2. Data Preprocessing
 Refer to the datapreprocess module to preprocess the raw single-cell multi-omics data, identify highly variable genes, and obtain the WNN (Weighted Nearest Neighbor) representation.

## 3. Preparation of Input Data
 Refer to the weight and motif modules, the WNN matrix is utilized to perform weighted integration of RNA and ATAC data, while also preparing motif data for subsequent transcription factor binding site analysis.

```
pbmc <- readRDS("pbmc_seurat_object.rds")
wnn_matrix <- readRDS("wnn_matrix.rds")
wnn_matrix_normalized <- wnn_matrix / rowSums(wnn_matrix) # Normalize each row (so that the sum of each row equals 1)

rna <- pbmc@assays$SCT@data
weighted_rna <- rna %*% wnn_matrix_normalized # RNA weighted average

atac <- pbmc@assays[["ATAC"]]@data # ATAC weighted average
weighted_atac <- atac %*% wnn_matrix_normalized

pwm_list <- chromVARmotifs::human_pwms_v2
pbmc_motif <- AddMotifs(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm_list)
```
## 4. CRE-gene cis-regulatory network inference
 Input processed seurat type data and highly variable genes, use runGenePeakLink function to obtain significant peak-gene relationships:

```
pbmc <- readRDS("pbmc_motif.rds")

file_path <- "highly_variable_genes.txt"
highly_variable_genes <- readLines(file_path) 

RNA <- GetAssayData(pbmc, assay = "SCT", layer = "data") 
RNA <- t(RNA)
peak_WNN <-readRDS("weighted_atac.rds")
peak_WNN <- t(peak_WNN)
peak_gene <- run_peak_gene(pbmc, RNA, peak_WNN, highly_variable_genes, lambda2 = 0.25, max_distance = 5e+05)
```
## 5. TF-gene regulatory network inference
  Input the seurat data and the peak-gene links obtained by filtering, and use the TFGeneGRN function to obtain the TF-gene relationship:

```
aen_lists <- make_aen_lists(peak_gene)
gene_names <- names(aen_lists)
RNA_WNN <- readRDS("weighted_rna.rds")
RNA_WNN <- RNA_WNN[gene_names, ]
RNA_WNN <- t(RNA_WNN)
tf_gene <- run_tf_gene(pbmc, RNA, peak_gene, RNA_WNN, gene_names, lambda2 = 0.5,layer = "data")
```