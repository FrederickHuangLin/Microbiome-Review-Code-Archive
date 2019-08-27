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
simparams = simparams%>%mutate(obs_seed = abn_seed+1)
simparams = simparams%>%separate(col = n_samp, into = c("n_samp_grp1", "n_samp_grp2"), sep = "_")
simparams = simparams%>%arrange(n_taxa, n_samp_grp1, prop_diff, abn_seed, obs_seed)
simparams_list = apply(simparams, 1, paste0, collapse = "_")

simparamslabels = c("n_taxa", "n_samp_grp1", "n_samp_grp2", "prop_diff", "abn_seed", "obs_seed")

source("scripts/ancom.R")
source("scripts/ancom_bc.R")
source("scripts/data_generation.R")

library(doParallel)
library(foreach)

detectCores()
myCluster = makeCluster(10, type = "FORK")
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
  pre_process = feature_table_pre_process(obs_abn, meta_data, 
                                          sample_var = "Sample_ID", group_var = "group",
                                          zero_cut = 0.90, lib_cut = 1000, neg_lb = FALSE)
  struc_zero = pre_process$structure_zeros
  num_struc_zero = apply(struc_zero, 1, sum)
  feature_table = pre_process$feature_table
  s0 = rownames(feature_table)[which(num_struc_zero == 0)]
  s1 = rownames(feature_table)[which(num_struc_zero == 1)]
  
  # Run ANCOM
  # Format for ANCOM: rows = subjects, cols=taxa
  feature_table_sub = t(feature_table[s0, ])
  
  res_W = ANCOM.main(OTUdat = feature_table_sub, Vardat = meta_data, 
                     main.var = "group", sig = 0.05, prev.cut = 1.01)
  
  res = data.frame(otu_names = c(s0, s1), W_stat = Inf, detected_0.9 = TRUE,
                   detected_0.8 = TRUE, detected_0.7 = TRUE, detected_0.6 = TRUE)
  res[match(as.character(res_W$otu.names), res$otu_names), ] = res_W
  res$diff_ind = test_dat$diff_taxa[c(s0, s1)]
  
  # FDR
  FDR = ifelse(sum(res$detected_0.7, na.rm = T) == 0, 0, 
               sum(ifelse(res$diff_ind == 0 & res$detected_0.7, 1, 0), na.rm = T)/
                 sum(res$detected_0.7, na.rm = T))
  
  # Power
  power = sum(ifelse(res$diff_ind != 0 & res$detected_0.7, 1, 0), na.rm = T)/
    sum(res$diff_ind != 0, na.rm = T)
  
  c(FDR, power)
}
end_time = Sys.time()
end_time - start_time

stopCluster(myCluster)
simlist = data.frame(simlist)
write_csv(simlist, "fdr_power_ancom.csv")