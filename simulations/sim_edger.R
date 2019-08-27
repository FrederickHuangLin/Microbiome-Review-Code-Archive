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

source("scripts/data_generation.R")

library(doParallel)
library(foreach)
library(edgeR)

start_time = Sys.time()
simlist = foreach(i = simparams_list, .combine = 'cbind') %do% {
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
  # Prepare data for edgeR
  groupdata = factor(rep(c(1, 2), c(n_samp_grp1, n_samp_grp2)))
  countdata = test_dat$obs_abn; meta_data = data.frame(group = groupdata)
  
  zero_threshold = 0.90
  taxa_info_ind = apply(countdata, 1, function(x) sum(x == 0)/(n_samp_grp1+n_samp_grp2))
  countdata = countdata[which(taxa_info_ind<zero_threshold), ] + 1L
  
  d = DGEList(counts = countdata, group = groupdata)
  d = calcNormFactors(d)
  design_mat = model.matrix(~ 0 + d$samples$group)
  colnames(design_mat) = levels(d$samples$group)
  
  d = estimateDisp(d, design_mat)
  fit = glmQLFit(d, design_mat)
  qlf = glmQLFTest(fit, contrast = c(1, -1))
  out = data.frame(taxa = rownames(topTags(qlf, n = nrow(countdata))$table), 
                   FDR = topTags(qlf, n = nrow(countdata))$table$FDR)
  out = out[match(rownames(countdata), as.character(out$taxa)), ]
  out$FDR[is.na(out$FDR)] = 1
  
  res = data.frame(diff_test = ifelse(out$FDR < 0.05, 1, 0), 
                   diff_ind = test_dat$diff_taxa[which(taxa_info_ind < zero_threshold)])
  
  # FDR
  FDR = ifelse(sum(res$diff_test == 1, na_rm = T) == 0, 0, 
               sum(ifelse(res$diff_ind == 0&res$diff_test == 1, 1, 0), na_rm = T)/
                 sum(res$diff_test == 1, na_rm = T))
  
  # Power
  power = sum(ifelse(res$diff_ind!= 0&res$diff_test == 1, 1, 0), na_rm = T)/
    sum(res$diff_ind!= 0, na_rm = T)
  
  c(FDR, power)
}
end_time = Sys.time()
end_time - start_time

simlist = data.frame(simlist)
write_csv(simlist, "fdr_power_edger.csv")