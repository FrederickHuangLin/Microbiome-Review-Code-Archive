library(tidyverse)
library(phyloseq)

abn_tab_gen = function(n_taxa, n_samp_grp1, n_samp_grp2, prop_diff, abn_seed, obs_seed, out_prop){
  # Total number of samples
  n_samp = n_samp_grp1 + n_samp_grp2
  
  set.seed(abn_seed) # This seed is used to control whether you would like to have the same population
  zero_prop = 0.2
  low_prop = 0.6 * (1 - zero_prop) # Proportion of low abundance 
  med_prop = 0.3 * (1 - zero_prop) # Proportion of medium abundance
  hi_prop = 0.1 * (1 - zero_prop)  # Proportion of high abundance
  # Indices for taxa abundance 
  index = sample(c(0, 1, 2, 3), n_taxa, replace = T, prob = c(zero_prop, low_prop, med_prop, hi_prop)) 
  
  # Poisson parameters for group 1
  lambda1 = rep(NA, n_taxa)
  lambda1[which(index==0)] = 0
  lambda1[which(index==1)] = rgamma(length(which(index==1)), shape=50, rate=1)
  lambda1[which(index==2)] = rgamma(length(which(index==2)), shape=200, rate=1)
  lambda1[which(index==3)] = rgamma(length(which(index==3)), shape=10000, rate=1)
  # Poisson parameters for group 2
  lambda2 = lambda1
  lambda2[which(index==0)] = sample(c(50, 200, 10000), sum(index==0), replace = T, prob = c(0.6, 0.3, 0.1))
  
  # Differentially abundant taxa
  diff_ind = rep(0, n_taxa)
  # Group1 is larger than group2
  diff1_ind = sample((1:n_taxa)[index!=0], floor(sum(index!=0)*prop_diff), replace=FALSE)
  diff_ind[diff1_ind] = 1
  # Group2 is higher than group1
  wt = runif(1, 0, 1)
  diff2_ind = sample(diff1_ind, wt*length(diff1_ind), replace=FALSE)
  diff_ind[diff2_ind] = 2
  # Structural zeros
  diff_ind[which(index==0)] = -1
  
  # Effect size
  effect_size = rep(1, n_taxa)
  effect_size[diff1_ind] = runif(length(diff1_ind), 1, 10)
  effect_size[diff2_ind] = runif(length(diff2_ind), 0.1, 1)
  names(effect_size) = paste0("taxon", seq(n_taxa))
  
  # Abundance template
  temp_grp1 = round(lambda1*effect_size)
  temp_grp2 = round(lambda2)
  temp_dat = data.frame(temp_grp1, temp_grp2, effect_size)
  rownames(temp_dat) = paste0("taxon", seq(n_taxa))
  
  # Abundance table in the ecosystem
  abn_mat = matrix(0, ncol = n_samp, nrow = n_taxa)
  for(i in 1:n_taxa){
    abn_mat[i, ] = c(rpois(n_samp_grp1, temp_grp1[i]), rpois(n_samp_grp2, temp_grp2[i]))
  }
  
  # Outliers
  out_ind = rep(0, n_taxa); out_ind[sample(seq(n_taxa), out_prop*n_taxa, replace  =  F)] = 1
  names(out_ind) = paste0("taxon", seq(n_taxa))
  abn_mat[which(out_ind == 1), sample(seq(n_samp), out_prop*n_samp, replace  =  F)] = 0
  
  abn_total = colSums(abn_mat)
  names(abn_total) = paste0("sub", seq(n_samp))
  
  # Number of taxa that are sampled for each subject
  depth = 1/sample(c(runif(n_samp, 10, 50), runif(n_samp, 100, 500)), n_samp, replace  =  T)
  obs_total = round(max(abn_total)*depth)
  names(obs_total) = paste0("sub", seq(n_samp))
  
  # Specimen abundance
  set.seed(obs_seed)
  obs_mat = 1:n_samp %>% map_dfc(function(i) 
    phyloseq:::rarefaction_subsample(x = abn_mat[, i], sample.size = obs_total[i]))
  
  # Prepare output data sets
  abn_dat = data.frame(abn_mat, row.names  =  NULL)
  rownames(abn_dat) = paste0("taxon", seq(n_taxa))
  colnames(abn_dat) = paste0("sub", seq(n_samp))
  
  obs_dat = data.frame(obs_mat, row.names  =  NULL)
  rownames(obs_dat) = paste0("taxon", seq(n_taxa))
  colnames(obs_dat) = paste0("sub", seq(n_samp))
  
  grp_ind = c(rep(1, n_samp_grp1), rep(2, n_samp_grp2))
  names(grp_ind) = paste0("sub", seq(n_samp))
  
  names(diff_ind) = paste0("taxon", seq(n_taxa))
  
  samp_frac = obs_total/abn_total
  names(samp_frac) = paste0("sub", seq(n_samp))
  
  test_data = list(temp_dat, abn_dat, obs_dat, effect_size, grp_ind, 
                   diff_ind, out_ind, samp_frac, abn_total, obs_total)
  names(test_data) = c("template", "pop_abn", "obs_abn", "effect_size", "grp", 
                       "diff_taxa", "outlier", "samp_frac", "abn_total", "obs_total")
  return(test_data)
}