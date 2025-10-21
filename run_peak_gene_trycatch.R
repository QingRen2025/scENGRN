run_peak_gene <- function(seurat,
                          pseudocell_mat,
                          pseudocell_mat_peaks,
                          gene_list,
                          species = "human",  # Add species as an argument
                          lambda2 = 0.5,
                          ci_cutoff = 1.96,
                          pval_cutoff = NULL,
                          set_nonsig_to_zero = TRUE,
                          max_distance = 5e+05,
                          train_fraction = 0.8,
                          num_bootstraps = 1000,
                          bootstrap = TRUE,
                          num_threads = 24,
                          globals_maxsize = NULL,
                          verbose = TRUE,
                          bs_seed = 147258,
                          multi_seed = 258369, ...) {
  
  require(gcdnet)
  require(future)
  require(future.apply)
  require(SeuratWrappers)
  require(Seurat)
  require(Signac)
  
  if (!is.null(gene_list)) {
    if (class(gene_list) == "list") {
      gene_list <- unlist(gene_list)
    }		
  } else {
    stop("Need a list of genes to model")
  }
  
  # Determine the gene coordinates and peak matrix based on the species
  if (species == "human") {
    # Code for human species (first function)
    gene_coords <- Signac:::CollapseToLongestTranscript(Annotation(seurat))
    gene_list <- gene_list[which(gene_list %in% gene_coords$gene_name)]
    
    peak_matrix <- Signac:::DistanceToTSS(peaks = StringToGRanges(colnames(pseudocell_mat_peaks)),
                                          genes = gene_coords[which(gene_coords$gene_name %in% gene_list)],
                                          distance = max_distance)
  } else if (species == "mouse") {
    # Code for mouse species (second function)
    annotations <- Annotation(seurat)
    annotations$gene_name <- str_to_title(annotations$gene_name)
    Annotation(seurat) <- annotations
    gene_coords <- Signac:::CollapseToLongestTranscript(Annotation(seurat))
    gene_list <- gene_list[which(gene_list %in% gene_coords$gene_name)]
    
    peak_matrix <- Signac:::DistanceToTSS(peaks = StringToGRanges(colnames(pseudocell_mat_peaks)),
                                          genes = gene_coords[which(gene_coords$gene_name %in% gene_list)],
                                          distance = max_distance)
  } else {
    stop("Species must be either 'human' or 'mouse'")
  }
  
  peak_matrix_list <- lapply(gene_list, function(x) {
    peak_col <- peak_matrix[, x]
    return(names(peak_col)[which(peak_col > 0)]) 
  })
  names(peak_matrix_list) <- gene_list
  
  # bootstrap setup
  num_bootstraps = num_bootstraps + 1 # last seed added for cross validation step
  if (!bootstrap) num_bootstraps = 2
  if (!(is.null(bs_seed))) set.seed(bs_seed)
  bootstrap_sequence_input_lists <- make_bootstrap_sequences(num_bootstraps, gene_list)
  
  # run aen
  peak_input_list <- lapply(gene_list, function(x) list(x, 
                                                        as.data.frame(pseudocell_mat_peaks[, peak_matrix_list[[x]]]), 
                                                        bootstrap_sequence_input_lists[[x]]))
  names(peak_input_list) <- gene_list
  
  # removing mitochondrial genes (no peaks to model)
  peak_input_list <- peak_input_list[which(unlist(lapply(peak_input_list, function(x) dim(x[[2]])[2])) > 1)]
  if (verbose && length(peak_input_list) != length(gene_list)) {
    print(paste0("Omitted (mitochondrial) genes with no genomic peaks: ",
                 gene_list[which(!(gene_list %in% names(peak_input_list)))]))
  } 
  
  if (num_threads > 1) {
    plan(multisession, workers = num_threads)
    if (!is.null(globals_maxsize)) options(future.globals.maxSize = globals_maxsize)
  }
  
  start <- Sys.time()
  if (num_threads == 1) {
    aen_results <- lapply(peak_input_list, run_peaks_for_gene,
                          pseudocell_matrix = pseudocell_mat,
                          train_fraction = train_fraction,
                          lambda2 = lambda2,
                          ci_cutoff = ci_cutoff,
                          pval_cutoff = pval_cutoff,
                          set_nonsig_to_zero = set_nonsig_to_zero, ...)
  } else if (!is.null(multi_seed)) {
    aen_results <- future_lapply(peak_input_list, run_peaks_for_gene,
                                 pseudocell_matrix = pseudocell_mat,
                                 train_fraction = train_fraction,
                                 lambda2 = lambda2,
                                 ci_cutoff = ci_cutoff,
                                 pval_cutoff = pval_cutoff,
                                 set_nonsig_to_zero = set_nonsig_to_zero,
                                 future.seed = multi_seed, ...)
  } else {
    aen_results <- future_lapply(peak_input_list, run_peaks_for_gene,
                                 pseudocell_matrix = pseudocell_mat,
                                 train_fraction = train_fraction,
                                 lambda2 = lambda2,
                                 ci_cutoff = ci_cutoff,
                                 pval_cutoff = pval_cutoff,
                                 set_nonsig_to_zero = set_nonsig_to_zero, ...)
  }
  names(aen_results) <- names(peak_input_list)
  end <- Sys.time()
  if (verbose) print(paste0("AEN completed in ", end - start))
  plan(sequential)
  
  # return object
  return(aen_results)
}

run_peaks_for_gene <- function(input_list,
                               pseudocell_matrix,
                               train_fraction = 0.8,
                               set_seed = TRUE,
                               lambda = "lambda.min",
                               lambda2 = 0.5,
                               ci_cutoff = 1.96,
                               pval_cutoff = NULL,
                               set_nonsig_to_zero = TRUE, ...) {
  require(gcdnet)
  
  target_name <- input_list[[1]]
  pseudocell_factor_expressions <- input_list[[2]]
  num_pseudo_cells <- dim(pseudocell_factor_expressions)[1]
  
  bootstrap_seeds <- input_list[[3]]
  if (set_seed && !is.null(bootstrap_seeds)) { 
    set.seed(bootstrap_seeds[length(bootstrap_seeds)])
  }
  
  train_rows <- sample(nrow(pseudocell_matrix), train_fraction * nrow(pseudocell_matrix))
  x.train <- as.matrix(pseudocell_factor_expressions[train_rows, ])
  x.test <- as.matrix(pseudocell_factor_expressions[-train_rows, ])
  
  
  y.train <- as.numeric(pseudocell_matrix[train_rows, target_name])
  y.test <- as.numeric(pseudocell_matrix[-train_rows, target_name])
  
  tryCatch({
    # 使用普通的弹性网络，而不是自适应弹性网络
    aen.fit <- cv.gcdnet(x.train, y.train, method = "ls", standardize = TRUE, 
                         intercept = TRUE, pred.loss = "loss", lambda2 = lambda2, 
                         ...)
    aen.predicted <- predict(aen.fit, s = lambda, newx = x.test)
    MSE <- 1/num_pseudo_cells * sum((y.test - aen.predicted)^2)
    TSS <- sum((y.test - mean(y.test))^2)
    r2 <- 1 - num_pseudo_cells * MSE/TSS
    num_bootstraps <- max(1, length(bootstrap_seeds) - 1)
    coef_df <- matrix(nrow = length(rownames(coef(aen.fit))), 
                      ncol = num_bootstraps)
    rownames(coef_df) <- rownames(coef(aen.fit))
    for (i in 1:num_bootstraps) {
      if (set_seed && !is.null(bootstrap_seeds)) 
        set.seed(bootstrap_seeds[i])
      bs_rows <- sample(nrow(pseudocell_matrix), nrow(pseudocell_matrix), 
                        replace = TRUE)
      x <- as.matrix(pseudocell_factor_expressions[bs_rows, 
      ])
      y <- as.matrix(pseudocell_matrix[bs_rows, target_name])
      lambda_to_use = ifelse(lambda == "lambda.min", aen.fit$lambda.min, 
                             aen.fit$lambda.1se)
      bs.aen <- gcdnet(x, y, method = "ls", standardize = TRUE, 
                       intercept = TRUE, lambda2 = lambda2, lambda = lambda_to_use, 
                       ...)
      coef_df[, i] <- coef(bs.aen)[rownames(coef_df), ]
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
    results <- list(target_name, r2, aen.fit, t(coefs))
    names(results) <- c("gene", "r2", "aen.fit", "coefs")
    return(results)
    
  }, error = function(e) {
    return(list(gene = target_name, r2 = NA, aen.fit = NULL, coefs = NULL))
  })
}