# Load dependencies
library(nloptr)
library(dplyr)

# Data Pre-Processing
feature_table_pre_process = function(feature_table, meta_data, sample_var, group_var, 
                                     zero_cut, lib_cut, neg_lb){
  feature_table = data.frame(feature_table, check.names = FALSE)
  meta_data = data.frame(meta_data, check.names = FALSE)
  # Drop unused levels
  meta_data[] = lapply(meta_data, function(x) if(is.factor(x)) factor(x) else x)
  
  sample_ID = colnames(feature_table)
  meta_data = meta_data[match(sample_ID, meta_data[, sample_var]), ]
  
  # 1. Identify outliers within each taxon
  group = factor(meta_data[, group_var])
  group_name = levels(group)
  grp_ind_origin = lapply(1:nlevels(group), function(i) which(group == group_name[i]))
  n_grp_origin = length(grp_ind_origin)
  n_samp_grp_origin = sapply(grp_ind_origin, length)
  feature_table = feature_table[, unlist(grp_ind_origin)]
  
  z = log(feature_table+1)
  f = z; f[f == 0] = NA; f = colMeans(f, na.rm = T)
  f_mean = unlist(tapply(f, rep(1:n_grp_origin, n_samp_grp_origin), mean))
  e = f-rep(f_mean, n_samp_grp_origin)
  y = t(t(z)-e)
  
  outlier_check = function(x){
    mu1 = quantile(x, 0.25); mu2 = quantile(x, 0.75)
    sigma1 = quantile(x, 0.75)-quantile(x, 0.25); sigma2 = sigma1
    pi = 0.75
    n = length(x)
    epsilon = 100
    tol = 1e-5
    score = pi*dnorm(x, mean = mu1, sd = sigma1)/((1-pi)*dnorm(x, mean = mu2, sd = sigma2))
    while (epsilon>tol) {
      grp1.ind = score >= 1
      mu1.new = mean(x[grp1.ind]); mu2_new = mean(x[!grp1.ind])
      sigma1.new = sd(x[grp1.ind]); if(is.na(sigma1.new)) sigma1.new = 0
      sigma2_new = sd(x[!grp1.ind]); if(is.na(sigma2_new)) sigma2_new = 0
      pi_new = sum(grp1.ind)/n
      
      para = c(mu1.new, mu2_new, sigma1.new, sigma2_new, pi_new)
      if(any(is.na(para))) break
      
      score = pi_new*dnorm(x, mean = mu1.new, sd = sigma1.new)/
        ((1-pi_new)*dnorm(x, mean = mu2_new, sd = sigma2_new))
      
      epsilon = sqrt((mu1-mu1.new)^2+(mu2-mu2_new)^2+
                     (sigma1-sigma1.new)^2+(sigma2-sigma2_new)^2+(pi-pi_new)^2)
      mu1 = mu1.new; mu2 = mu2_new; sigma1 = sigma1.new; sigma2 = sigma2_new; pi = pi_new
    }
    
    if(mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2){
      if(pi > 0.85){
        out_ind = (!grp1.ind)
      }else if(pi < 0.15){
        out_ind = grp1.ind
      }else{
        out_ind = rep(FALSE, n)
      }
    }else{
      out_ind = rep(FALSE, n)
    }
    return(out_ind)
  }
  feature_table_out = t(apply(y, 1, function(i)
    unlist(tapply(i, rep(1:n_grp_origin, n_samp_grp_origin), function(j) outlier_check(j)))))
  feature_table[feature_table_out] = NA
  
  # 2. Discard taxa with zeros >=  pre_cut
  taxa_zero_prop = apply(feature_table, 1, function(x) sum(x == 0, na.rm = T)/length(x[!is.na(x)]))
  filter_taxa = which(taxa_zero_prop >= zero_cut)
  if(length(filter_taxa)>0){
    feature_table = feature_table[-filter_taxa, ]
  }
  
  # 3. Discard samples with library size < lib_cut
  library_size = colSums(feature_table, na.rm = T)
  sample_ID = colnames(feature_table)
  meta_data = meta_data[match(sample_ID, meta_data[, sample_var]), ]
  
  if(any(library_size < lib_cut)){
    filter_subject = which(library_size < lib_cut)
    feature_table = feature_table[, - filter_subject]
    meta_data = meta_data[- filter_subject, ]
  }
  
  # 4. Re-order the OTU table
  group = factor(meta_data[, group_var])
  group_name = levels(group)
  grp_ind = lapply(1:nlevels(group), function(i) which(group == group_name[i]))
  n_grp = length(grp_ind)
  n_samp_grp = sapply(grp_ind, length)
  
  n_taxa = nrow(feature_table)
  taxa_id = rownames(feature_table)
  n_samp = ncol(feature_table)
  
  # 5. Identify taxa with structure zeros
  present_table = as.matrix(feature_table)
  present_table[is.na(present_table)] = 0
  present_table[present_table != 0] = 1
  
  p_hat_mat = t(apply(present_table, 1, function(x)
    unlist(tapply(x, rep(1:n_grp, n_samp_grp), function(y) mean(y, na.rm = T)))))
  sample_size = t(apply(feature_table, 1, function(x)
    unlist(tapply(x, rep(1:n_grp, n_samp_grp), function(y) length(y[!is.na(y)])))))
  p_hat_lo_mat = p_hat_mat-1.96*sqrt(p_hat_mat*(1-p_hat_mat)/sample_size)
  colnames(p_hat_mat) = levels(group)
  colnames(p_hat_lo_mat) = levels(group)
  
  struc_zero = matrix(0, nrow = n_taxa, ncol = n_grp)
  struc_zero[p_hat_mat == 0] = 1
  # Whether we need to classify a taxon into structural zero by its negative lower bound?
  if(neg_lb) struc_zero[p_hat_lo_mat <= 0] = 1
  rownames(struc_zero) = taxa_id
  colnames(struc_zero) = paste0("Structural Zero in ", levels(group))
  
  # 6_ Return results
  res = list(feature_table = feature_table, library_size = library_size, 
             group_name = group_name, group_ind = grp_ind, structure_zeros = struc_zero)
  return(res)
}

# ANCOM-BC main function
ANCOM_BC = function(feature_table, grp_name, grp_ind, struc_zero, adj_method, 
                  tol_EM, max_iterNum, alpha){
  n_taxa_origin = nrow(feature_table)
  taxa_id_origin = rownames(feature_table)
  n_samp = ncol(feature_table)
  sample_id = colnames(feature_table)
  
  n_grp = length(grp_ind)
  n_samp_grp = sapply(grp_ind, length)
  
  ### 0. Separate out taxa with no structural zeros
  info_taxa_pos = which(apply(struc_zero, 1, function(x) all(x == 0)))
  O = feature_table[info_taxa_pos, ]
  n_taxa = nrow(O)
  taxa_id = rownames(O)
  n_samp = ncol(O)
  y = log(O+1)
  
  ### 1. Estimate sampling fractions and mean abundances by iteratively least squares
  d = rep(0, n_samp)
  mu = t(apply(y, 1, function(i) tapply(i, rep(1:n_grp, n_samp_grp), function(j)
    mean(j, na.rm = T))))
  iterNum = 0
  epsilon = 100
  while (epsilon > tol_EM & iterNum < max_iterNum) {
    # Updating mu
    mu_new = t(apply(t(t(y) - d), 1, function(i) tapply(i, rep(1:n_grp, n_samp_grp), function(j)
      mean(j, na.rm = T))))
    
    # Updating d
    d_new = colMeans(y - mu_new[, rep(1:ncol(mu_new), times = n_samp_grp)], na.rm = T)
    
    # Iteration
    epsilon = sqrt(sum((mu_new - mu)^2) + sum((d_new-d)^2))
    iterNum = iterNum+1
    
    mu = mu_new
    d = d_new
  }
  
  mu_var_raw = (y-t(t(mu[, rep(1:ncol(mu), times = n_samp_grp)])+d))^2
  mu_var = t(apply(mu_var_raw, 1, function(x) tapply(x, rep(1:n_grp, n_samp_grp), function(y)
    mean(y, na.rm = T))))
  sample_size = t(apply(y, 1, function(x)
    unlist(tapply(x, rep(1:n_grp, n_samp_grp), function(y) length(y[!is.na(y)])))))
  mu_var = mu_var/sample_size
  
  ### 2_ Estimate the bias of sampling fractions by E-M algorithm
  Delta = mu[, 1]-mu[, 2]
  Delta_var_est = rowSums(mu_var)
  
  ## 2_1 Initials
  pi1.0 = 0.125
  pi2_0 = 0.75
  pi3_0 = 0.125
  delta_0 = mean(Delta[Delta >= quantile(Delta, 0.25, na.rm = T)&
                       Delta <= quantile(Delta, 0.75, na.rm = T)], na.rm = T)
  d1.0 = mean(Delta[Delta < quantile(Delta, 0.125, na.rm = T)], na.rm = T)
  d2_0 = mean(Delta[Delta > quantile(Delta, 0.875, na.rm = T)], na.rm = T)
  psi1.sq_0 = var(Delta[Delta < quantile(Delta, 0.125, na.rm = T)], na.rm = T)
  if(is.na(psi1.sq_0)|psi1.sq_0 == 0) psi1.sq_0 = 1
  psi2_sq_0 = var(Delta[Delta>quantile(Delta, 0.875, na.rm = T)], na.rm = T)
  if(is.na(psi2_sq_0)|psi2_sq_0 == 0) psi2_sq_0 = 1
  
  ## 2_2 Apply E-M algorithm
  # 2_21 Store all paras in vectors/matrices
  pi1.vec = c(pi1.0); pi2_vec = c(pi2_0); pi3_vec = c(pi3_0)
  delta_vec = c(delta_0); d1.vec = c(d1.0); d2_vec = c(d2_0)
  psi1.sq_vec = c(psi1.sq_0); psi2_sq_vec = c(psi2_sq_0)
  
  # 2_22 E-M iteration
  iterNum = 0
  epsilon = 100
  sigmai_sq = Delta_var_est
  while (epsilon > tol_EM & iterNum < max_iterNum) {
    # print(iterNum)
    ## Current value of paras
    pi1 = pi1.vec[length(pi1.vec)]; pi2 = pi2_vec[length(pi2_vec)]; pi3 = pi3_vec[length(pi3_vec)]
    delta = delta_vec[length(delta_vec)]; 
    d1 = d1.vec[length(d1.vec)]; d2 = d2_vec[length(d2_vec)]
    psi1.sq = psi1.sq_vec[length(psi1.sq_vec)]; psi2_sq = psi2_sq_vec[length(psi2_sq_vec)]
    
    ## E-step
    pdf1 = sapply(seq(n_taxa), function(i) dnorm(Delta[i], delta+d1, sqrt(sigmai_sq[i]+psi1.sq)))
    pdf2 = sapply(seq(n_taxa), function(i) dnorm(Delta[i], delta, sqrt(sigmai_sq[i])))
    pdf3 = sapply(seq(n_taxa), function(i) dnorm(Delta[i], delta+d2, sqrt(sigmai_sq[i]+psi2_sq)))
    r1i = pi1 * pdf1/(pi1 * pdf1+pi2 * pdf2+pi3 * pdf3); r1i[is.na(r1i)] = 0
    r2i = pi2 * pdf2/(pi1 * pdf1+pi2 * pdf2+pi3 * pdf3); r2i[is.na(r2i)] = 0
    r3i = pi3 * pdf3/(pi1 * pdf1+pi2 * pdf2+pi3 * pdf3); r3i[is.na(r3i)] = 0
    
    ## M-step
    pi1.new = mean(r1i, na.rm = T); pi2_new = mean(r2i, na.rm = T); pi3_new = mean(r3i, na.rm = T)
    delta_new = sum(r1i * (Delta - d1) / (sigmai_sq + psi1.sq) + r2i * Delta / sigmai_sq + 
                      r3i * (Delta - d2) / (sigmai_sq + psi2_sq), na.rm = T) / 
      sum(r1i / (sigmai_sq + psi1.sq) + r2i / sigmai_sq + r3i / (sigmai_sq + psi2_sq), na.rm = T)
    d1.new = min(sum(r1i * (Delta-delta) / (sigmai_sq + psi1.sq), na.rm = T) / 
                 sum(r1i / (sigmai_sq + psi1.sq), na.rm = T), 0)
    if(is.nan(d1.new)) d1.new = 0
    d2_new = max(sum(r3i * (Delta-delta) / (sigmai_sq + psi2_sq), na.rm = T) / 
                 sum(r3i / (sigmai_sq + psi2_sq), na.rm = T), 0)
    if(is.nan(d2_new)) d2_new = 0
    
    # Nelder-Mead simplex algorithm for psi1.sq, psi2_sq, and sigmai_sq
    obj_psi1.sq = function(x){
      log_pdf = log(sapply(seq(n_taxa), function(i) dnorm(Delta[i], delta + d1, sqrt(sigmai_sq[i] + x))))
      log_pdf[is.infinite(log_pdf)] = 0
      - sum(r1i*log_pdf, na.rm = T)
    }
    psi1.sq_new = neldermead(x0 = psi1.sq, fn = obj_psi1.sq, lower = 0)$par
    
    obj_psi2_sq = function(x){
      log_pdf = log(sapply(seq(n_taxa), function(i) dnorm(Delta[i], delta + d2, sqrt(sigmai_sq[i] + x))))
      log_pdf[is.infinite(log_pdf)] = 0
      - sum(r3i*log_pdf, na.rm = T)
    }
    psi2_sq_new = neldermead(x0 = psi2_sq, fn = obj_psi2_sq, lower = 0)$par
    
    ## Merge to the paras vectors/matrices
    pi1.vec = c(pi1.vec, pi1.new); pi2_vec = c(pi2_vec, pi2_new); pi3_vec = c(pi3_vec, pi3_new)
    delta_vec = c(delta_vec, delta_new)
    d1.vec = c(d1.vec, d1.new); d2_vec = c(d2_vec, d2_new)
    psi1.sq_vec = c(psi1.sq_vec, psi1.sq_new); psi2_sq_vec = c(psi2_sq_vec, psi2_sq_new)
    
    ## Calculate the new epsilon
    epsilon = sqrt((pi1.new - pi1)^2 + (pi2_new - pi2)^2 + (pi3_new - pi3)^2 + 
                     (delta_new - delta)^2 + (d1.new - d1)^2 + 
                     (d2_new - d2)^2 + (psi1.sq_new - psi1.sq)^2 + (psi2_sq_new - psi2_sq)^2)
    iterNum = iterNum + 1
  }
  bias_est = delta_vec[length(delta_vec)]
  
  ### 3_ Test results
  ## 3_1 Results for taxa with non-structural zeros
  W_numerator = mu[, 1] - mu[, 2] - bias_est
  W_denominator = mu_var[, 1] + mu_var[, 2]
  
  W = W_numerator/sqrt(W_denominator)
  p_val = sapply(W, function(x) 2*pnorm(abs(x), mean = 0, sd = 1, lower.tail = F))
  q_val = p.adjust(p_val, method = adj_method)
  q_val[is.na(q_val)] = 1
  
  res_nonstrc_zero = data.frame(W_numerator, se = sqrt(W_denominator), W, p_val, q_val)
  
  ## 3_2 Results for taxa with structural zeros
  if(length(info_taxa_pos) < n_taxa_origin){
    O_strc_zero = feature_table[-info_taxa_pos, ]
    ind_strc_zero = struc_zero[-info_taxa_pos, rep(1:n_grp, times = n_samp_grp)]
    O_strc_zero_adj = O_strc_zero*(1 - ind_strc_zero)
    y_strc_zero = log(O_strc_zero_adj+1)
    d_strc_zero = t(t(1 - ind_strc_zero)*d)
    y_strc_zero_adj = y_strc_zero - d_strc_zero
    mu_strc_zero = t(apply(y_strc_zero_adj, 1, function(i) 
      tapply(i, rep(1:n_grp, n_samp_grp), function(j) mean(j, na.rm = T))))
    # Make it the relative mean difference (with related to the smallest value)
    mu_strc_zero_adj = mu_strc_zero
    mu_strc_zero_adj[mu_strc_zero_adj == 0] = NA
    mu_strc_zero_adj = t(t(mu_strc_zero_adj) + abs(apply(mu_strc_zero, 2, min)))
    mu_strc_zero_adj[is.na(mu_strc_zero_adj)] = 0
    
    res_strc_zero = data.frame(W_numerator = mu_strc_zero_adj[, 1] - mu_strc_zero_adj[, 2], 
                               se = 0, W = Inf, p_val = 0, q_val = 0)
  }else{
    res_strc_zero = NA
  }
  
  ## 3_3 Combine results together
  res = data.frame(W_numerator = Inf, se = 0, W = Inf, 
                   p_val = rep(0, n_taxa_origin), q_val = rep(0, n_taxa_origin))
  res[info_taxa_pos, ] = res_nonstrc_zero
  res[-info_taxa_pos, ] = res_strc_zero
  res = mutate(res, diff_abn = ifelse(q_val < alpha, TRUE, FALSE))
  
  colnames(res) = c(paste0("mean_difference (", grp_name[1], " - ", grp_name[2], ")"), 
                    "se", paste0("effect_size (", grp_name[1], " - ", grp_name[2], ")"), 
                    "p_val", "q_val", "diff_abn")
  out = list(feature_table = feature_table, res = res, d = d, mu = mu, bias_est = bias_est)
  return(out)
}
