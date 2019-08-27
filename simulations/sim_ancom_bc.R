library(tidyverse)

# The number of taxa, sampling depth, and sample size
n_taxa = 1000; n_samp = c("20_30", "50_50")

# The proportion of differentially abundant taxa
prop_diff = c(0.05, 0.15, 0.25)

# Set seeds
iterNum = 100
abn_seed = seq(iterNum)

# Define the simulation parameters combinations
simparams = expand.grid(n_taxa, n_samp, prop_diff, abn_seed)
colnames(simparams) = c("n_taxa", "n_samp", "prop_diff", "abn_seed")
simparams = simparams%>%mutate(obs_seed = abn_seed + 1)
simparams = simparams%>%separate(col = n_samp, into = c("n_samp_grp1", "n_samp_grp2"), sep = "_")
simparams = simparams%>%arrange(n_taxa, n_samp_grp1, prop_diff, abn_seed, obs_seed)
simparams_list = apply(simparams, 1, paste0, collapse = "_")

simparamslabels = c("n_taxa", "n_samp_grp1", "n_samp_grp2", "prop_diff", "abn_seed", "obs_seed")

source("scripts/ancom_bc.R")
source("scripts/data_generation.R")

library(doParallel)
library(foreach)

detectCores()
myCluster = makeCluster(2, type = "FORK")
registerDoParallel(myCluster)

start_time = Sys.time()
simlist = foreach(i = simparams_list, .combine = 'cbind') %dopar% {
  # i = simparams_list[1]
  print(i)
  params = strsplit(i, "_")[[1]]
  names(params) = simparamslabels
  
  # Paras for data generation
  n_taxa = as.numeric(params["n_taxa"])
  n_samp_grp1 = as.numeric(params["n_samp_grp1"])
  n_samp_grp2 = as.numeric(params["n_samp_grp2"])
  prop_diff = as.numeric(params["prop_diff"])
  abn_seed = as.numeric(params["abn_seed"])
  obs_seed = as.numeric(params["obs_seed"])
  
  # Data generation
  test_dat = abn_tab_gen(n_taxa, n_samp_grp1, n_samp_grp2, prop_diff, abn_seed, obs_seed,
                         out_prop = 0.05)
  obs_abn = test_dat$obs_abn
  meta_data = cbind(Sample_ID = paste0("sub", seq(n_samp_grp1+n_samp_grp2)), 
                    group = rep(c(1, 2), c(n_samp_grp1, n_samp_grp2)))
  
  # Pre-processing
  feature_table = obs_abn; sample_var = "Sample_ID"; group_var = "group"; 
  zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
  pre_process = feature_table_pre_process(feature_table, meta_data, sample_var, 
                                          group_var, zero_cut, lib_cut, neg_lb)
  feature_table = pre_process$feature_table
  group_name = pre_process$group_name
  group_ind = pre_process$group_ind
  struc_zero = pre_process$structure_zeros
  
  # Paras for ANCOM-BC
  grp_name = group_name; grp_ind = group_ind; adj_method = "bonferroni"
  tol_EM = 1e-5; max_iterNum = 100; alpha = 0.05
  
  # Run ANCOM-BC
  suppressWarnings(out <- try(ANCOM_BC(feature_table, grp_name, grp_ind, struc_zero,
                                       adj_method, tol_EM, max_iterNum, alpha), 
                              silent = TRUE))
  if (inherits(out, "try-error")) {
    FDR = NA; power = NA
  }else{
    res = cbind(out$res, diff_ind = test_dat$diff_taxa[rownames(out$feature_table)])
    
    # FDR
    FDR = ifelse(sum(res$diff_abn, na.rm = T) == 0, 0, 
                 sum(ifelse(res$diff_ind == 0&res$diff_abn, 1, 0), na.rm = T)/
                   sum(res$diff_abn, na.rm = T))
    
    # Power
    power = sum(ifelse(res$diff_ind!= 0&res$diff_abn, 1, 0), na.rm = T)/
      sum(res$diff_ind!= 0, na.rm = T)
  }
  c(FDR, power)
}
end_time = Sys.time()
end_time - start_time

stopCluster(myCluster)
simlist = data.frame(simlist)
write_csv(simlist, "fdr_power_ancom_bc.csv")