run_tf_gene <- function(seurat,
                        RNA,
                        en_results_peaks = peak_results,
                        spatial_lag_RNA = spatial_lag_RNA,
                        gene_list = NULL,
                        lambda2 = 0.5,
                        ci_cutoff = 1.96,
                        pval_cutoff = NULL,
                        set_nonsig_to_zero = TRUE,
                        promoter_only = FALSE,
                        promoter_threshold = 2000,
                        train_fraction = 0.8,
                        num_bootstraps = 100,
                        bootstrap = TRUE,
                        num_threads = 24,
                        globals_maxsize = NULL,
                        verbose = TRUE,
                        bs_seed = 147258,
                        peak_assay = "ATAC",
                        layer = "data",
                        multi_seed = 258369, ...) {
  require(Seurat)
  require(Signac)
  require(SeuratWrappers)
  require(BiocGenerics)
  require(tidyverse)
  require(gcdnet)
  require(future)
  require(future.apply)
  
  if (!is.null(gene_list) && !promoter_only) {
    if (class(gene_list) == "list") {
      gene_list <- unlist(gene_list)
    }		
    modeled_genes <- gene_list
  } else if (!promoter_only) {
    modeled_genes <- names(en_results_peaks)
  } else {
    stop("Need a list of genes to model")
  }
  
  peak_tf_key <- seurat@assays$ATAC@motifs@data
  colnames(peak_tf_key) <- unlist(seurat@assays$ATAC@motifs@motif.names)
  
  regulator_tf_names <- colnames(peak_tf_key)
  regulator_tf_names <- regulator_tf_names[which(regulator_tf_names %in% colnames(RNA))]
  
  peak_tf_key <- peak_tf_key[, regulator_tf_names]
  peak_tf_key_tibble <- as_tibble(BiocGenerics::t(peak_tf_key), rownames = "rowname")
  
  # prepare initial TF-gene lists
  if (DefaultAssay(seurat) != peak_assay) DefaultAssay(seurat) <- peak_assay
  gene_coords <- Signac:::CollapseToLongestTranscript(Annotation(seurat))
  all_peaks <- rownames(GetAssayData(seurat, assay = peak_assay, layer = layer))
  prom_matrix <- upstreamDistanceToTSS(peaks = StringToGRanges(all_peaks),
                                       genes = gene_coords[which(gene_coords$gene_name %in% modeled_genes)],
                                       distance = promoter_threshold)
  
  prom_input_list <- lapply(modeled_genes, function(x) {
    peak_col <- prom_matrix[, x]
    return(names(peak_col)[which(peak_col > 0)])
  })
  names(prom_input_list) <- modeled_genes
  
  #prom_input_list <- prom_input_list[which(names(prom_input_list) %in% names(peak_input_list))] # get rid of same 3 MT genes
  
  if (promoter_only) {
    # linked_prom_peak_list <- lapply(prom_input_list, function(x) colnames(x[[2]]))
    tf_aggregates <- lapply(prom_input_list, peaks_to_tfs, key_tibble = peak_tf_key_tibble)
    names(tf_aggregates) <- modeled_genes
    tf_aggregates <- tf_aggregates[which(unlist(lapply(tf_aggregates, length)) > 1)]
    genes_without_annotated_promoters <- modeled_genes[which(!(modeled_genes %in% names(tf_aggregates)))]
    if (verbose) print(paste0("Genes without any candidate TFs--most likely no peaks in promoter region based on threshold and provided genome annotation: ", 
                              genes_without_annotated_promoters))
  } else {
    en_lists <- lapply(en_results_peaks, function(x) {
      l <- x[[4]]
      l <- as.data.frame(l)[which(l[, "coef_if_kept"] != 0), ]
      l <- rownames(l)[which(rownames(l) != "(Intercept)")]
    })
    names(en_lists) <- modeled_genes
    linked_peak_list <- lapply(modeled_genes, function(x) {
      union(en_lists[[x]], prom_input_list[[x]])
    })
    names(linked_peak_list) <- modeled_genes
    
    tf_aggregates <- lapply(linked_peak_list, peaks_to_tfs, key_tibble = peak_tf_key_tibble)
    names(tf_aggregates) <- modeled_genes
    tf_aggregates <- tf_aggregates[which(unlist(lapply(tf_aggregates, length)) > 1)]
    genes_without_peaks <- modeled_genes[which(!(modeled_genes %in% names(tf_aggregates)))]
    if (verbose) print(paste0("Genes without any candidate TFs--most likely no linked peaks or promoter region peaks: ",
                              genes_without_peaks))
  }
  
  # bootstrap setup
  num_bootstraps = num_bootstraps + 1 # last seed added for cross validation step
  if (!bootstrap) num_bootstraps = 2
  if (!(is.null(bs_seed))) set.seed(bs_seed)
  bootstrap_sequence_input_lists <- make_bootstrap_sequences(num_bootstraps, gene_list)
  
  print(num_bootstraps)
  print(bootstrap_sequence_input_lists[[1]][[1]])
  
  # run en
  tf_aggregates <- tf_aggregates[which(unlist(lapply(tf_aggregates, length)) > 1)] 
  input_list <- lapply(names(tf_aggregates), function(x) list(x, 
                                                              tf_aggregates[[x]],
                                                              bootstrap_sequence_input_lists[[x]]))
  names(input_list) <- names(tf_aggregates)
  
  if (num_threads > 1) {
    plan(multisession, workers = num_threads)
    if (!is.null(globals_maxsize)) options(future.globals.maxSize = globals_maxsize)
  }
  start <- Sys.time()
  if (num_threads == 1) {
    en_results <- lapply(input_list, run_tfs_for_gene,
                          regulator_names = regulator_tf_names,
                          RNA_matrix = RNA,
                          spatial_lag_RNA = spatial_lag_RNA,
                          lambda2 = lambda2,
                          train_fraction = train_fraction,
                          ci_cutoff = ci_cutoff,
                          pval_cutoff = pval_cutoff,
                          set_nonsig_to_zero = set_nonsig_to_zero, ...)
  } else if (!is.null(multi_seed)) {
    en_results <- future_lapply(input_list, run_tfs_for_gene,
                                 regulator_names = regulator_tf_names,
                                 RNA_matrix = RNA,
                                 spatial_lag_RNA = spatial_lag_RNA,
                                 lambda2 = lambda2,
                                 train_fraction = train_fraction,
                                 ci_cutoff = ci_cutoff,
                                 pval_cutoff = pval_cutoff,
                                 set_nonsig_to_zero = set_nonsig_to_zero,
                                 future.seed = multi_seed, ...)
  } else {
    en_results <- future_lapply(input_list, run_tfs_for_gene,
                                 regulator_names = regulator_tf_names,
                                 RNA_matrix = RNA,
                                 spatial_lag_RNA = spatial_lag_RNA,
                                 lambda2 = lambda2,
                                 train_fraction = train_fraction,
                                 ci_cutoff = ci_cutoff,
                                 pval_cutoff = pval_cutoff,
                                 set_nonsig_to_zero = set_nonsig_to_zero, ...)
  }
  names(en_results) <- names(input_list)
  end <- Sys.time()
  if (verbose) print(paste0("en completed in ", end - start))
  plan(sequential)
  
  return(en_results)
}

#################################################################
upstreamDistanceToTSS <- function(peaks,
                                  genes,
                                  distance = 200000,
                                  sep = c("-", "-")
) {
  require(Matrix)
  require(GenomicRanges)
  
  tss <- resize(x = genes, width = 1, fix = 'start')
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = 0
    )
  )
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  )
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  )
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- genes.extended$gene_name
  return(hit_matrix)
}

run_tfs_for_gene <- function (input_tf_list, RNA_matrix, spatial_lag_RNA, regulator_names, 
                              train_fraction = 0.8, set_seed = TRUE, lambda = "lambda.min", 
                              lambda2 = 0.5, ci_cutoff = 1.96, pval_cutoff = NULL, 
                              set_nonsig_to_zero = TRUE, ...) 
{
  require(gcdnet)
  target_name <- input_tf_list[[1]]
  tf_agg <- input_tf_list[[2]]
  names(tf_agg) <- tf_agg
  
  if (is.null(spatial_lag_RNA)) {
    stop("spatial_lag_RNA must be provided")
  }
  # Get the spatial lag column for the target gene (ensure column names match)
  spatial_lag_col <- spatial_lag_RNA[, target_name, drop = FALSE]
  colnames(spatial_lag_col) <- paste0(target_name, "_spatial_lag")  # Add suffix to avoid naming conflicts
  
  # Ensure rownames are aligned (critical step)
  spatial_lag_col <- spatial_lag_col[rownames(RNA_matrix), , drop = FALSE]
  
  if (is.null(regulator_names)) {
    tf_agg <- tf_agg[which(tf_agg %in% colnames(RNA_matrix))]
  }
  else {
    tf_agg <- tf_agg[which(tf_agg %in% regulator_names)]
  }
  
  # Combine spatial lag values when creating the feature matrix
  pseudocell_factor_expressions <- cbind(
    RNA_matrix[, tf_agg],  # Original TF expression matrix
    spatial_lag_col        # Newly added spatial lag column
  )
  colnames(pseudocell_factor_expressions) <- c(tf_agg, colnames(spatial_lag_col))
  
  num_pseudo_cells <- dim(RNA_matrix)[1]
  bootstrap_seeds <- input_tf_list[[3]]
  if (set_seed && !is.null(bootstrap_seeds)) {
    set.seed(bootstrap_seeds[length(bootstrap_seeds)])
  }
  if (target_name %in% tf_agg) {
    pseudocell_factor_expressions <- pseudocell_factor_expressions[, 
                                                                   (names(pseudocell_factor_expressions) != target_name)]
  }
  train_rows = sample(nrow(RNA_matrix), train_fraction * 
                        nrow(RNA_matrix))
  x.train = as.matrix(pseudocell_factor_expressions[train_rows, 
  ])
  x.test = as.matrix(pseudocell_factor_expressions[-train_rows, 
  ])
  y.train = as.numeric(RNA_matrix[train_rows, target_name])
  y.test = as.numeric(RNA_matrix[-train_rows, target_name])
  
  tryCatch({
    en.fit <- cv.gcdnet(x.train, y.train, method = "ls", standardize = TRUE, 
                         intercept = TRUE, pred.loss = "loss", lambda2 = lambda2, 
                         ...)
    en.predicted <- predict(en.fit, s = lambda, newx = x.test)
    MSE <- 1/num_pseudo_cells * sum((y.test - en.predicted)^2)
    TSS <- sum((y.test - mean(y.test))^2)
    r2 <- 1 - num_pseudo_cells * MSE/TSS
    num_bootstraps <- max(1, length(bootstrap_seeds) - 1)
    coef_df <- matrix(nrow = length(rownames(coef(en.fit))), 
                      ncol = num_bootstraps)
    rownames(coef_df) <- rownames(coef(en.fit))
    for (i in 1:num_bootstraps) {
      if (set_seed && !is.null(bootstrap_seeds)) 
        set.seed(bootstrap_seeds[i])
      bs_rows <- sample(nrow(RNA_matrix), nrow(RNA_matrix), 
                        replace = TRUE)
      x <- as.matrix(pseudocell_factor_expressions[bs_rows, 
      ])
      y <- as.matrix(RNA_matrix[bs_rows, target_name])
      lambda_to_use = ifelse(lambda == "lambda.min", en.fit$lambda.min, 
                             en.fit$lambda.1se)
      bs.en <- gcdnet(x, y, method = "ls", standardize = TRUE, 
                       intercept = TRUE, lambda2 = lambda2, lambda = lambda_to_use, 
                       ...)
      coef_df[, i] <- coef(bs.en)[rownames(coef_df), ]
    }
    if (set_nonsig_to_zero && is.null(ci_cutoff)) {
      ci_cutoff <- qnorm(1 - pval_cutoff/2, 0, 1)
    }
    coefs <- apply(coef_df, 1, function(x) {
      coef_mean <- mean(x)
      if (length(x) > 1 && num_bootstraps > 1) {
        standard_error <- sqrt(1/(length(x) - 1) * sum((x - 
                                                          coef_mean)^2))
        se_cutoff <- ci_cutoff * standard_error
      }
      else {
        standard_error <- NA
        se_cutoff <- 0
      }
      if (set_nonsig_to_zero) {
        c(ifelse(abs(coef_mean) > se_cutoff, coef_mean, 0), 
          coef_mean, standard_error)
      }
      else {
        c(coef_mean, coef_mean, standard_error)
      }
    }, simplify = TRUE)
    if (is.null(dim(coefs))) {
      names(coefs) <- c("coef_if_kept", "coef_mean", "se")
    }
    else {
      rownames(coefs) <- c("coef_if_kept", "coef_mean", "se")
    }
    results <- list(target_name, r2, en.fit, t(coefs))
    names(results) <- c("gene", "r2", "en.fit", "coefs")
    return(results)
    
  }, error = function(e) {
    return(list(gene = target_name, r2 = NA, en.fit = NULL, coefs = NULL))
  })
}