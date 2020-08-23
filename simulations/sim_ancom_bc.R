library(tidyverse)
library(doParallel)
library(foreach)

detectCores()
myCluster = makeCluster(4, type = "FORK")
registerDoParallel(myCluster)

source("data_generation.R")
source("ancom_bc.R")

n_taxa = 200; n_samp = 60
x = data.frame(group = paste0("G", rep(1:2, each = n_samp/2))); type = "none"; group = "group"
prop_diff = c(0.05, 0.15, 0.25); zero_prop = 0.2; depth = "small"
meta_data = data.frame(sample_id = paste0("sample", seq(n_samp)), x)

# Set seeds
iterNum = 100
abn_seed = seq(iterNum)

# Define the simulation parameters
simparams = expand.grid(prop_diff, abn_seed)
colnames(simparams) = c("prop_diff", "abn_seed")
simparams = simparams %>% mutate(obs_seed = abn_seed + 1) %>% arrange(prop_diff, abn_seed, obs_seed)
simparams_list = apply(simparams, 1, paste0, collapse = "_")
simparamslabels = c("prop_diff", "abn_seed", "obs_seed")

simlist = foreach(i = simparams_list, .combine = 'cbind', 
                  .packages = c("MASS", "nloptr", "CVXR")) %dopar% {
  # i = simparams_list[[1]]
  print(i)
  params = strsplit(i, "_")[[1]]
  names(params) = simparamslabels
  
  # Paras for data generation
  prop_diff = as.numeric(params["prop_diff"])
  abn_seed = as.numeric(params["abn_seed"])
  obs_seed = as.numeric(params["obs_seed"])
  
  # Data generation
  test_dat = abn_tab_gen(n_taxa, n_samp, x, type, group, prop_diff,
                         abn_seed, obs_seed, zero_prop, depth)
  obs_abn = test_dat$obs_abn
  
  # Run ANCOM-BC
  feature_table = obs_abn; sample_id = "sample_id"; p_adj_method = "holm"
  zero_cut = 0.90; lib_cut = 0; tol = 1e-5; max_iter = 100; conserve = FALSE
  alpha = 0.05; per_num = 1000; adj_formula = "group"; struc_zero = TRUE; neg_lb = FALSE
  global = FALSE; direct = FALSE; dunnett = FALSE; pattern = NULL
  
  suppressWarnings(out <- try(ANCOM_BC(feature_table, meta_data, sample_id, adj_formula,
                                       p_adj_method, zero_cut, lib_cut, struc_zero, neg_lb, group, 
                                       tol, max_iter, conserve, alpha, per_num, 
                                       global, direct, dunnett, pattern), silent = TRUE))
  
  if (inherits(out, "try-error")) {
    FDR = NA; power = NA
  }else{
    res_test = out$res$diff_abn[, 2] * 1
    res_true = test_dat$diff_ind * 1
    res_true[test_dat$zero_ind] = 1
    res_true = res_true[rownames(out$feature_table)]
    TP = sum(res_test[res_true == 1] == 1, na.rm = T)
    FP = sum(res_test[res_true == 0] == 1, na.rm = T)
    FN = sum(res_test[res_true == 1] == 0, na.rm = T)
    FDR = FP/(TP + FP); power = TP/(TP + FN)
  }
  c(FDR, power)
}

stopCluster(myCluster)
write_csv(data.frame(simlist), "sim_fdr_power_ancom_bc.csv")

