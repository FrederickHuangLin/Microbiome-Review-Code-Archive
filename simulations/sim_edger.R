library(tidyverse)
library(doParallel)
library(foreach)
library(edgeR)

source("data_generation.R")

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
                  .packages = c("edgeR")) %do% {
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
  
  # Prepare data
  countdata = test_dat$obs_abn
  rownames(meta_data) = meta_data$sample_id
  
  zero_threshold = 0.90
  taxa_info_ind = apply(countdata, 1, function(x) sum(x == 0)/n_samp)
  feature_table = round(countdata[which(taxa_info_ind < zero_threshold), ]) + 1L
  
  d = DGEList(counts = feature_table, group = meta_data$group)
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
  
  res_test = ifelse(out$FDR < 0.05, 1, 0)
  res_true = test_dat$diff_ind * 1
  res_true[test_dat$zero_ind] = 1
  res_true = res_true[rownames(feature_table)]
  TP = sum(res_test[res_true == 1] == 1, na.rm = T)
  FP = sum(res_test[res_true == 0] == 1, na.rm = T)
  FN = sum(res_test[res_true == 1] == 0, na.rm = T)
  FDR = FP/(TP + FP); power = TP/(TP + FN)
  c(FDR, power)
}

write_csv(data.frame(simlist), "sim_fdr_power_edger.csv")