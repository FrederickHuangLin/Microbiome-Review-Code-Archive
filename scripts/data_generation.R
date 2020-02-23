library(phyloseq)

abn_tab_gen = function(n_taxa, n_samp, x, type = NULL, group = NULL, prop_diff, 
                       abn_seed, obs_seed, zero_prop = 0.2, depth){
  set.seed(abn_seed) # The seed to control the population
  low_prop = 0.6 # Proportion of low abundance 
  med_prop = 0.3 # Proportion of medium abundance
  hi_prop = 0.1  # Proportion of high abundance
  # Indices for absolute abundance 
  index = sample(1:3, n_taxa, replace = T, prob = c(low_prop, med_prop, hi_prop)) 
  
  # Covariates of interest
  covariates = colnames(x)
  x = data.frame(x, check.names = FALSE)
  x[] = lapply(x, function(i) if(is.numeric(i)) scale(i, center = TRUE, scale = TRUE) else i)
  x_formula = formula(paste0("~", paste(covariates, collapse = " + ")))
  x_expand = model.matrix(x_formula, data = x)
  
  # Coefficients in log scale
  # Overall mean absolute abundance
  b0 = rep(NA, n_taxa)
  b0[which(index == 1)] = log(50)
  b0[which(index == 2)] = log(200)
  b0[which(index == 3)] = log(10000)
  # The matrix of coefficients
  b = list(b0 = b0)
  # Taxa of differntially abundant
  diff_pos = sample(1:n_taxa, prop_diff * n_taxa)
  # Node for the specified order
  node = matrix(NA, nrow = n_taxa, ncol = length(covariates))
  colnames(node) = covariates
  rownames(node) = paste0("taxon", 1:n_taxa)
  for (k in 1:length(covariates)) {
    b_k = matrix(0, nrow = n_taxa, 
                 ncol = ncol(x_expand[, grepl(covariates[k], colnames(x_expand)), drop = FALSE]))
    if (type[k] == "none") {
      # No order
      wt = runif(1, 0, 1)
      diff_pos1 = sample(diff_pos, wt * length(diff_pos), replace = FALSE)
      diff_pos2 = setdiff(diff_pos, diff_pos1)
      b_k[diff_pos1, ] = log(runif(length(diff_pos1) * ncol(b_k), 1, 10))
      b_k[diff_pos2, ] = log(runif(length(diff_pos2) * ncol(b_k), 0.1, 1))
      node[, k] = rep(0, n_taxa)
    } else if (type[k] == "simple") {
      # Simple order (monotonic increasing)
      b_k[diff_pos, ] = log(runif(length(diff_pos) * ncol(b_k), 1, 10))
      # Sort the coefficients only when the number of levels > 2
      if (ncol(b_k) > 1) b_k = t(apply(b_k, 1, sort))
      # Nodes
      node_k = rep(0, n_taxa)
      node_k[diff_pos] = ncol(b_k)
      node[, k] = node_k
    } else if (type[k] == "tree") {
      # Tree order (>= the reference group)
      b_k[diff_pos, ] = log(runif(length(diff_pos) * ncol(b_k), 1, 10))
      # Nodes
      node_k = rep(0, n_taxa)
      node_k[diff_pos] = 1
      node[, k] = node_k
    } else {
      # Umbrella order (first increasing then decreasing)
      b_k[diff_pos, ] = log(runif(length(diff_pos) * ncol(b_k), 1, 10))
      # Define the function for generating umbrella order
      umbrella = function(x, node) {
        x_max = max(x); max_pos = which.max(x)
        x_node = x[node]
        x[max_pos] = x_node
        x[node] = x_max
        x[1:(node - 1)] = sort(x[1:(node - 1)])
        x[(node + 1):length(x)] = sort(x[(node + 1):length(x)], decreasing = TRUE)
        return(x)
      }
      max_pos = ceiling(ncol(b_k)/2)
      if (ncol(b_k) > 1) {
        # Apply the function only when the number of levels > 2
        b_k = t(apply(b_k, 1, function(x) umbrella(x, max_pos)))
      }
      # Nodes
      node_k = rep(0, n_taxa)
      node_k[diff_pos] = max_pos
      node[, k] = node_k
    }
    b = c(b, list(b_k))
  }
  diff_ind = rep(FALSE, n_taxa)
  diff_ind[diff_pos] = TRUE
  names(diff_ind) = paste0("taxon", 1:n_taxa)
  b = Reduce('cbind', b)
  colnames(b) = colnames(x_expand)
  rownames(b) = paste0("taxon", 1:n_taxa)
  
  # Error term
  std = sqrt(1/exp(b0))
  names(std) = paste0("taxon", 1:n_taxa)
  
  # Log absoulte abundance in the ecosystem
  log_abn_mat = matrix(NA, ncol = n_samp, nrow = n_taxa)
  for (i in 1:n_taxa) {
    log_abn_mat[i, ] = x_expand %*% b[i, ] + rnorm(n_samp, mean = 0, sd = std[i])
  }
  colnames(log_abn_mat) = paste0("sample", 1:n_samp)
  rownames(log_abn_mat) = paste0("taxon", 1:n_taxa)
  
  # Absolute abundance in the ecosystem
  abn_mat = ceiling(exp(log_abn_mat))
  abn_mat = abn_mat - 1
  colnames(abn_mat) = paste0("sample", 1:n_samp)
  rownames(abn_mat) = paste0("taxon", 1:n_taxa)
  
  # Structural zeros
  if (!is.null(group)) {
    zero_ind = sort(sample(setdiff(1:n_taxa, diff_pos), zero_prop * n_taxa))
    # Set the absolute abundance of the reference group to be 0
    # For simplicity, only set structural zeros on the first level of the first group variable
    ref_group = colnames(x_expand)[grepl(group[1], colnames(x_expand))][1]
    abn_mat[zero_ind, which(x_expand[, ref_group] == 1)] = 0
    b[zero_ind, ref_group] = - b[zero_ind, 1]
  }
  
  # Microbial load
  micro_load = colSums(abn_mat)
  names(micro_load) = paste0("sample", seq(n_samp))
  
  # Library size
  if(depth == "large"){
    depth = 1/runif(n_samp, 5, 10)
  }else{
    depth = 1/sample(c(runif(n_samp, 10, 50), runif(n_samp, 100, 500)), n_samp, replace = T)
  }
  lib_size = round(max(micro_load) * depth)
  names(lib_size) = paste0("sample", 1:n_samp)
  
  # Absolute abundance in the sample
  set.seed(obs_seed)
  obs_list = lapply(1:n_samp, function(i) 
    phyloseq:::rarefaction_subsample(x = abn_mat[, i], sample.size = lib_size[i]))
  obs_mat = Reduce('cbind', obs_list)
  colnames(obs_mat) = paste0("sample", 1:n_samp)
  rownames(obs_mat) = paste0("taxon", 1:n_taxa)
  
  # Sampling fraction
  c = lib_size/micro_load
  names(c) = paste0("sample", 1:n_samp)
  
  # Outputs
  test_data = list(b, x, std, diff_ind, zero_ind, node, abn_mat, obs_mat, c, micro_load, lib_size)
  names(test_data) = c("coeff", "covariates", "std", "diff_ind", "zero_ind", "node",
                       "pop_abn", "obs_abn", "samp_frac", "micro_load", "lib_size")
  return(test_data)
}