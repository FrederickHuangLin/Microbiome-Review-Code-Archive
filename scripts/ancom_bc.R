library(dplyr)
library(MASS)
library(nloptr)
library(CVXR)

# OTU table should be a matrix/data.frame with each feature in rows and sample in columns. 
# Metadata should be a matrix/data.frame containing the sample identifier. 

ANCOM_BC = function(feature_table, meta_data, sample_id, adj_formula, p_adj_method = "holm", 
                    zero_cut = 0.90, lib_cut = 1000, struc_zero = FALSE, neg_lb = FALSE, group = NULL, 
                    tol = 1e-5, max_iter = 100, conserve = FALSE, alpha = 0.05, per_num = 1000,
                    global = FALSE, direct = FALSE, dunnett = FALSE, pattern = NULL){
  # 0. Data preprocessing
  feature_table = data.frame(feature_table, check.names = FALSE)
  meta_data = data.frame(meta_data, check.names = FALSE)
  # Drop unused levels
  meta_data[] = lapply(meta_data, function(x) if(is.factor(x)) factor(x) else x)
  # Match sample IDs between metadata and feature table
  id = intersect(meta_data[, sample_id], colnames(feature_table))
  feature_table = feature_table[, id]
  meta_data = meta_data[match(id, meta_data[, sample_id]), ]
  # Check group variables
  if (!is.null(group)) {
    n_level = apply(meta_data[, group, drop = FALSE], 2, function(x) length(unique(x)))
    if (any(n_level < 2)) stop("Group variables should have >=2 levels.")
  }
  
  ## 0.1 Discard taxa with zeros  >=  zero_cut
  zero_prop = apply(feature_table, 1, function(x) sum(x == 0, na.rm = T)/length(x[!is.na(x)]))
  taxa_del = which(zero_prop >= zero_cut)
  if(length(taxa_del) > 0){
    feature_table = feature_table[- taxa_del, ]
  }
  
  ## 0.2 Discard samples with library size < lib_cut
  lib_size = colSums(feature_table, na.rm = T)
  if(any(lib_size < lib_cut)){
    subj_del = which(lib_size < lib_cut)
    feature_table = feature_table[, - subj_del]
    meta_data = meta_data[- subj_del, ]
  }

  n_taxa = nrow(feature_table)
  taxa_id = rownames(feature_table)
  n_samp = ncol(feature_table)
  samp_id = colnames(feature_table)
  
  ## 0.3 Identify taxa with structural zeros
  if (struc_zero) {
    if (is.null(group)) stop("Please specify the group variables for identifying structural zeros.")
    zero_ind = vector(mode = "list", length = length(group))
    names(zero_ind) = group
    for (k in 1:length(group)) {
      group_k = factor(meta_data[, group[k]])
      present_table = as.matrix(feature_table)
      present_table[is.na(present_table)] = 0
      present_table[present_table != 0] = 1
      
      p_hat = t(apply(present_table, 1, function(x)
        unlist(tapply(x, group_k, function(y) mean(y, na.rm = T)))))
      samp_size = t(apply(feature_table, 1, function(x)
        unlist(tapply(x, group_k, function(y) length(y[!is.na(y)])))))
      p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)
      
      zero_ind_k = (p_hat == 0) * 1
      # Whether we need to classify a taxon into structural zero by its negative lower bound?
      if (neg_lb) zero_ind_k[p_hat_lo <= 0] = 1
      
      colnames(zero_ind_k) = paste0("structural_zero (", 
                                        group[k], " = ",
                                        colnames(zero_ind_k), ")")
      zero_ind[[k]] = zero_ind_k
    }
  }else{
    zero_ind = NULL
  }
  
  ## 0.4 OTU table transformation: 
  # Add pseudocount (1) and take logarithm.
  y = log(feature_table + 1)
  options(na.action = "na.pass") # Keep NA's in rows
  x = model.matrix(formula(paste0("~", adj_formula)), data = meta_data)
  options(na.action = "na.omit") # Change it back
  covariates = colnames(x)
  n_covariates = length(covariates)
  
  # 1. Estimation of parameters 
  ## 1.1 Sampling fractions
  d = rep(0, n_samp)
  tformula = formula(paste0("y ~ ", adj_formula))
  fits = lapply(1:n_taxa, function(i) {
    df = data.frame(y = unlist(y[i, ]) - d, meta_data)
    return(lm(tformula, data = df))
  })
  ## 1.2 Regression coefficients
  beta = lapply(fits, function(i) {
    beta_i = rep(NA, length(covariates)) # prevent errors of missing values
    coef_i = coef(i)
    beta_i[match(names(coef_i), covariates)] = coef_i
    return(beta_i)
  })
  beta = Reduce('rbind', beta)
  
  ## 1.3 Iterative least square
  iterNum = 0
  epsilon = 100
  while (epsilon > tol & iterNum < max_iter) {
    # Updating beta
    fits = lapply(1:n_taxa, function(i) {
      df = data.frame(y = unlist(y[i, ]) - d, meta_data)
      return(lm(tformula, data = df))
    })
    beta_new = lapply(fits, function(i) {
      beta_i = rep(NA, length(covariates))
      coef_i = coef(i)
      beta_i[match(names(coef_i), covariates)] = coef_i
      return(beta_i)
    })
    beta_new = Reduce('rbind', beta_new)

    # Updating d
    y_hat = lapply(fits, function(i) {
      y_hat_i = rep(NA, n_samp)
      fit_i = fitted(i)
      y_hat_i[match(names(fit_i), samp_id)] = fit_i
      return(y_hat_i)

    })
    y_hat = Reduce('rbind', y_hat)
    d_new = colMeans(y - y_hat, na.rm = T)

    # Iteration
    epsilon = sqrt(sum((beta_new - beta)^2, na.rm = T) + sum((d_new - d)^2, na.rm = T))
    iterNum = iterNum + 1

    beta = beta_new
    d = d_new
  }
  colnames(beta) = covariates
  rownames(beta) = taxa_id
  names(d) = samp_id
  
  ## 1.4 Regression residuals
  y_hat = lapply(fits, function(i) {
    y_hat_i = rep(NA, n_samp)
    fit_i = fitted(i)
    y_hat_i[match(names(fit_i), samp_id)] = fit_i
    return(y_hat_i)
    
  })
  y_hat = Reduce('rbind', y_hat)
  colnames(y_hat) = samp_id
  rownames(y_hat) = taxa_id
  e = t(t(y - y_hat) - d)
  
  ## 1.5 Variance-covariance matrices of coefficients
  XTX_inv = ginv(t(x[complete.cases(x), ]) %*% x[complete.cases(x), ])
  var_cov_hat = vector(mode = "list", length = n_taxa) # Covariances
  var_hat = matrix(NA, nrow = n_taxa, ncol = n_covariates) # Variances
  for (i in 1:n_taxa) {
    sigma2_xxT = matrix(0, ncol = n_covariates, nrow = n_covariates)
    for (j in 1:n_samp) {
      sigma2_xxT_j = e[i, j]^2 * x[j, ] %*% t(x[j, ])
      sigma2_xxT_j[is.na(sigma2_xxT_j)] = 0
      sigma2_xxT = sigma2_xxT + sigma2_xxT_j
    }
    var_cov_hat[[i]] = XTX_inv %*% sigma2_xxT %*% XTX_inv
    rownames(var_cov_hat[[i]]) = covariates
    colnames(var_cov_hat[[i]]) = covariates
    var_hat[i, ] = diag(var_cov_hat[[i]])
  }
  names(var_cov_hat) = taxa_id
  colnames(var_hat) = covariates
  rownames(var_hat) = taxa_id
  
  # 2. Estimation of the between-group bias by E-M algorithm
  delta_em = rep(NA, ncol(beta) - 1)
  delta_wls = rep(NA, ncol(beta) - 1)
  var_delta = rep(NA, ncol(beta) - 1)
  for (i in 1:length(delta_em)) {
    # Ignore the intercept
    Delta = beta[, i + 1]; Delta = Delta[!is.na(Delta)]
    nu0 = var_hat[, i + 1]; nu0 = nu0[!is.na(nu0)]
    
    ## 2.1 Initials
    pi0_0 = 0.75
    pi1_0 = 0.125
    pi2_0 = 0.125
    delta_0 = mean(Delta[Delta >= quantile(Delta, 0.25, na.rm = T)&
                           Delta <= quantile(Delta, 0.75, na.rm = T)], na.rm = T)
    l1_0 = mean(Delta[Delta < quantile(Delta, 0.125, na.rm = T)], na.rm = T)
    l2_0 = mean(Delta[Delta > quantile(Delta, 0.875, na.rm = T)], na.rm = T)
    kappa1_0 = var(Delta[Delta < quantile(Delta, 0.125, na.rm = T)], na.rm = T)
    if(is.na(kappa1_0)|kappa1_0 == 0) kappa1_0 = 1
    kappa2_0 = var(Delta[Delta > quantile(Delta, 0.875, na.rm = T)], na.rm = T)
    if(is.na(kappa2_0)|kappa2_0 == 0) kappa2_0 = 1
    
    ## 2.2 Apply E-M algorithm
    ### 2.21 Store all paras in vectors/matrices
    pi0_vec = c(pi0_0); pi1_vec = c(pi1_0); pi2_vec = c(pi2_0)
    delta_vec = c(delta_0); l1_vec = c(l1_0); l2_vec = c(l2_0)
    kappa1_vec = c(kappa1_0); kappa2_vec = c(kappa2_0)
    
    ### 2.22 E-M iteration
    iterNum = 0
    epsilon = 100
    while (epsilon > tol & iterNum < max_iter) {
      # Current value of paras
      pi0 = pi0_vec[length(pi0_vec)]; pi1 = pi1_vec[length(pi1_vec)]; pi2 = pi2_vec[length(pi2_vec)]
      delta = delta_vec[length(delta_vec)]; 
      l1 = l1_vec[length(l1_vec)]; l2 = l2_vec[length(l2_vec)]
      kappa1 = kappa1_vec[length(kappa1_vec)]; kappa2 = kappa2_vec[length(kappa2_vec)]
      
      # E-step
      pdf0 = sapply(seq(n_taxa), function(i) dnorm(Delta[i], delta, sqrt(nu0[i])))
      pdf1 = sapply(seq(n_taxa), function(i) dnorm(Delta[i], delta + l1, sqrt(nu0[i] + kappa1)))
      pdf2 = sapply(seq(n_taxa), function(i) dnorm(Delta[i], delta + l2, sqrt(nu0[i] + kappa2)))
      r0i = pi0*pdf0/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2); r0i[is.na(r0i)] = 0
      r1i = pi1*pdf1/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2); r1i[is.na(r1i)] = 0
      r2i = pi2*pdf2/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2); r2i[is.na(r2i)] = 0
      
      # M-step
      pi0_new = mean(r0i, na.rm = T); pi1_new = mean(r1i, na.rm = T); pi2_new = mean(r2i, na.rm = T)
      delta_new = sum(r0i*Delta/nu0 + r1i*(Delta-l1)/(nu0+kappa1) + r2i*(Delta-l2)/(nu0+kappa2), na.rm = T)/
        sum(r0i/nu0 + r1i/(nu0+kappa1) + r2i/(nu0+kappa2), na.rm = T)
      l1_new = min(sum(r1i*(Delta-delta)/(nu0+kappa1), na.rm = T)/sum(r1i/(nu0+kappa1), na.rm = T), 0)
      l2_new = max(sum(r2i*(Delta-delta)/(nu0+kappa2), na.rm = T)/sum(r2i/(nu0+kappa2), na.rm = T), 0)
      
      # Nelder-Mead simplex algorithm for kappa1 and kappa2
      obj_kappa1 = function(x){
        log_pdf = log(sapply(seq(n_taxa), function(i) dnorm(Delta[i], delta+l1, sqrt(nu0[i]+x))))
        log_pdf[is.infinite(log_pdf)] = 0
        -sum(r1i*log_pdf, na.rm = T)
      }
      kappa1_new = neldermead(x0 = kappa1, fn = obj_kappa1, lower = 0)$par
      
      obj_kappa2 = function(x){
        log_pdf = log(sapply(seq(n_taxa), function(i) dnorm(Delta[i], delta+l2, sqrt(nu0[i]+x))))
        log_pdf[is.infinite(log_pdf)] = 0
        -sum(r2i*log_pdf, na.rm = T)
      }
      kappa2_new = neldermead(x0 = kappa2, fn = obj_kappa2, lower = 0)$par
      
      # Merge to the paras vectors/matrices
      pi0_vec = c(pi0_vec, pi0_new); pi1_vec = c(pi1_vec, pi1_new); pi2_vec = c(pi2_vec, pi2_new)
      delta_vec = c(delta_vec, delta_new)
      l1_vec = c(l1_vec, l1_new); l2_vec = c(l2_vec, l2_new)
      kappa1_vec = c(kappa1_vec, kappa1_new); kappa2_vec = c(kappa2_vec, kappa2_new)
      
      # Calculate the new epsilon
      epsilon = sqrt((pi0_new-pi0)^2 + (pi1_new-pi1)^2 + (pi2_new-pi2)^2 + (delta_new-delta)^2+
                       (l1_new-l1)^2 + (l2_new-l2)^2 + (kappa1_new-kappa1)^2 + (kappa2_new-kappa2)^2)
      iterNum = iterNum+1
    }
    ### 2.23 The EM estimator of bias
    delta_em[i] = delta_vec[length(delta_vec)]
    
    ### 2.24 The WLS estimator of bias
    # Cluster 0
    C0 = which(Delta >= quantile(Delta, pi1_new, na.rm = T) & Delta < quantile(Delta, 1 - pi2_new, na.rm = T))
    # Cluster 1
    C1 = which(Delta < quantile(Delta, pi1_new, na.rm = T))
    # Cluster 2
    C2 = which(Delta >= quantile(Delta, 1 - pi2_new, na.rm = T))
    # Numerator of the WLS estimator
    nu = nu0
    nu[C1] = nu[C1] + kappa1_new
    nu[C2] = nu[C2] + kappa2_new
    wls_deno = sum(1 / nu)
    # Denominator of the WLS estimator
    wls_nume = 1 / nu
    wls_nume[C0] = (wls_nume * Delta)[C0]
    wls_nume[C1] = (wls_nume * (Delta - l1_new))[C1]
    wls_nume[C2] = (wls_nume * (Delta - l2_new))[C2]   
    wls_nume = sum(wls_nume)
    
    delta_wls[i] = wls_nume / wls_deno
    
    ### 2.25 Estimate the variance of bias  
    var_delta[i] = 1 / wls_deno
    if (is.na(var_delta[i])) var_delta[i] = 0
  }
  
  # 3. Coefficients, standard error, and sampling fractions
  beta_hat = beta
  beta_hat[, -1] = t(t(beta_hat[, -1]) - delta_em)
  
  if (conserve) {
    # Account for the variance of delta_hat
    se_hat = sqrt(sweep(var_hat, 2, c(0, var_delta), "+") + 
                    2 * sqrt(sweep(var_hat, 2, c(0, var_delta), "*")))
  }else{
    se_hat = sqrt(var_hat)
  }
  
  d_hat = vector()
  for (i in 1:n_taxa) {
    d_hat_i = y[i, ] - x %*% beta_hat[i, ]
    d_hat = rbind(d_hat, d_hat_i)
  }
  d_hat = colMeans(d_hat, na.rm = TRUE)
  
  # 4. Test results
  ## 4.1 Model summary
  W = beta_hat/se_hat
  p = 2 * pnorm(abs(W), mean = 0, sd = 1, lower.tail = FALSE)
  q = apply(p, 2, function(x) p.adjust(x, method = p_adj_method))
  diff_abn = ifelse(q < alpha, TRUE, FALSE)
  res = list(beta = beta_hat, se = se_hat, W = W, 
             p_val = p, q_val = q, diff_abn = diff_abn)
  
  # Activate multi-group comparison only when there are group variables with >= 3 levels
  group_all = group
  zero_ind_all = zero_ind
  if (all(n_level < 3)) {
    global = FALSE; direct = FALSE; dunnett = FALSE; pattern = NULL
  } else if (any(n_level < 3)) {
    keep_ind = which(n_level >= 3)
    group = group_all[keep_ind]
    zero_ind = zero_ind_all[keep_ind]
    pattern = lapply(pattern, function(x) x[keep_ind])
  }
  
  ## 4.2 Global test
  # Define the function of global test
  global_test = function(para, p_adj_method){
    W_global = matrix(NA, nrow = n_taxa, ncol = length(para))
    p_global = matrix(NA, nrow = n_taxa, ncol = length(para))
    for (k in 1:length(para)) {
      k_ind = grepl(para[k], covariates)
      # Loop over the parameter of interest
      beta_hat_k = beta_hat[, k_ind]
      var_cov_hat_k = lapply(var_cov_hat, function(x) x = x[k_ind, k_ind])
      for (i in 1:n_taxa) {
        # Loop over taxa
        beta_hat_ik = beta_hat_k[i, ]
        var_cov_hat_ik = var_cov_hat_k[[i]]
        A = diag(x = 1, nrow = length(beta_hat_ik))
        W = t(A %*% beta_hat_ik) %*% ginv(A %*% var_cov_hat_ik %*% t(A)) %*% (A %*% beta_hat_ik)
        p = 2 * min(pchisq(W, df = length(beta_hat_ik), lower.tail = TRUE), 
                    pchisq(W, df = length(beta_hat_ik), lower.tail = FALSE))
        W_global[i, k] = W
        p_global[i, k] = p
      }
    }
    # Model summary
    rownames(W_global) = taxa_id; colnames(W_global) = para
    rownames(p_global) = taxa_id; colnames(p_global) = para
    q_global = apply(p_global, 2, function(x) p.adjust(x, method = p_adj_method))
    q_global[is.na(q_global)] = 1
    diff_global = ifelse(q_global < alpha, TRUE, FALSE)
    
    res_global = list(W = W_global, p_val = p_global, q_val = q_global, diff_abn = diff_global)
    return(res_global)
  }
  
  # Global test summary
  if (global) {
    if (is.null(group)) stop("Global test requires specifying the group variables.") 
    res_global = global_test(para = group, p_adj_method = p_adj_method)
  } else {
    res_global = NULL
  }
  
  ## 4.3 Directional test
  # Define the function of combination
  combn_fun = function(x, fun) {
    y = c(x, combn(x, 2, FUN = fun))
    combn_name = combn(names(x), 2)
    combn_name = paste(combn_name[2, ], combn_name[1, ], sep = " - ")
    names(y) = c(names(x), combn_name)
    return(y)
  }
  # Define the function for pairwise comparison
  pair_test = function(para) {
    beta_hat_pair = vector(mode = "list", length = length(para))
    names(beta_hat_pair) = para
    # Standard error for pairwise difference of parameters of interest
    var_hat_pair = vector(mode = "list", length = length(para))
    names(var_hat_pair) = para
    # Run the combination function to obtain pairwise comparison results
    for (k in 1:length(para)) {
      k_ind = grepl(para[k], covariates)
      beta_hat_k = beta_hat[, k_ind]
      var_hat_k = var_hat[, k_ind]
      beta_hat_pair_k = t(apply(beta_hat_k, 1, function(x) combn_fun(x, fun = diff)))
      var_hat_pair_k = t(apply(var_hat_k, 1, function(x) combn_fun(x, fun = sum)))
      beta_hat_pair[[k]] = beta_hat_pair_k
      var_hat_pair[[k]] = var_hat_pair_k
    }
    beta_hat_pair = Reduce('cbind', beta_hat_pair)
    var_hat_pair = Reduce('cbind', var_hat_pair)
    return(list(beta_hat_pair = beta_hat_pair, var_hat_pair = var_hat_pair))
  }
  
  # Directional test summary
  if (direct) {
    if (is.null(group)) stop("Directional test requires specifying the group variables.") 
    if (dunnett) {
      # For Dunnett's type of test: test against the reference group
      beta_hat_direct = beta_hat[, grepl(group, covariates), drop = FALSE]
      var_hat_direct = var_hat[, grepl(group, covariates), drop = FALSE]
    } else {
      # For general pairwise comparison
      # Coefficients for pairwise difference of parameters of interest
      res_direct = pair_test(para = group)
      beta_hat_direct = res_direct$beta_hat_pair 
      var_hat_direct = res_direct$var_hat_pair
    }
    W_direct = beta_hat_direct/sqrt(var_hat_direct)
    
    # mdFDR controlling procedure
    # Run the global test
    res_screen = global_test(para = group, p_adj_method = "BH")
    # The total number of null hypotheses rejected in the global test
    R = colSums(res_screen$diff_abn)
    # P-values for pairwise tests
    p_direct = 2 * pnorm(abs(W_direct), mean = 0, sd = 1, lower.tail = FALSE)
    # Only consider R significant taxa with regards to the global test
    n_mult = sapply(group, function(x) sum(grepl(x, colnames(p_direct))))
    screen_ind = res_screen$diff_abn[, rep(seq(length(group)), n_mult)]
    p_direct = p_direct * screen_ind
    p_direct[p_direct == 0] = 1
    # Adjust pairwise p-values at level of R/m*alpha
    q_direct = vector(mode = "list", length = length(group))
    for (k in 1:length(group)) {
      k_ind = grepl(group[k], colnames(p_direct))
      p_direct_k = p_direct[, k_ind, drop = FALSE]
      q_direct_k = t(apply(p_direct_k, 1, function(x) 
        p.adjust(x, method = p_adj_method, n = length(x) * n_taxa /R[k])))
      q_direct[[k]] = q_direct_k
    }
    q_direct = Reduce('cbind', q_direct)
    
    # Directional test summary
    diff_direct = ifelse(q_direct < alpha, TRUE, FALSE)
    se_hat_direct = sqrt(var_hat_direct)
    res_direct = list(beta = beta_hat_direct, se = se_hat_direct, W = W_direct, 
                      p_val = p_direct, q_val = q_direct, diff_abn = diff_direct)
    } else {
    res_direct = NULL
  }
  
  # 4.4 Test for patterns
  constrain_est = function(beta_hat, var_cov_hat, A, increase) {
    beta_opt = vector()
    for (i in 1:n_taxa) {
      beta_hat_i = beta_hat[i, ]
      var_cov_hat_i = var_cov_hat[[i]]
      # Estimated coefficients under constraint
      beta_opt_i = Variable(rows = length(beta_hat_i), cols = 1)
      obj = Minimize(matrix_frac(beta_opt_i - beta_hat_i, var_cov_hat_i))
      if (increase) {
        cons = A %*% beta_opt_i >= 0
      } else {
        cons = A %*% beta_opt_i <= 0
      }
      problem = Problem(objective = obj, constraints = list(cons))
      result = solve(problem)
      beta_opt_i = as.numeric(result$getValue(beta_opt_i))
      # Treat missing value as 0
      if (any(is.na(beta_opt_i))) beta_opt_i = rep(0, length(beta_hat_i))
      beta_opt = rbind(beta_opt, beta_opt_i)
    }
    rownames(beta_opt) = taxa_id; colnames(beta_opt) = colnames(beta_hat)
    return(beta_opt)
  }
  
  trend_test = function(para, pattern) {
    beta_opt = vector(mode = "list", length = length(para))
    W_opt = matrix(NA, nrow = n_taxa, ncol = length(para))
    p_opt = matrix(NA, nrow = n_taxa, ncol = length(para))
    q_opt = matrix(NA, nrow = n_taxa, ncol = length(para))
    node = matrix(NA, nrow = n_taxa, ncol = length(para))
    rownames(node) = taxa_id; colnames(node) = para
    rownames(W_opt) = taxa_id; colnames(W_opt) = para
    rownames(p_opt) = taxa_id; colnames(p_opt) = para
    rownames(q_opt) = taxa_id; colnames(q_opt) = para
    for (k in 1:length(para)) {
      k_ind = grepl(para[k], covariates)
      beta_hat_k = beta_hat[, k_ind]
      var_cov_hat_k = lapply(var_cov_hat, function(x) x = x[k_ind, k_ind])
      var_hat_k = var_hat[, k_ind]
      # Estimated coefficients under null
      beta_null_k = vector()
      for (i in 1:n_taxa) {
        var_cov_hat_ik = var_cov_hat_k[[i]]
        beta_null_k = rbind(beta_null_k, 
                            mvrnorm(n = per_num, 
                                    mu = rep(0, sum(k_ind)), 
                                    Sigma = var_cov_hat_ik))
      }
      if (!pattern[k] %in% c("simple", "tree", "umbrella")) {
        stop("Type of order should be either simple, tree, or umbrella.")
      }
      if (pattern[k] == "simple") {
        # Simple order
        A = diag(1, nrow = sum(k_ind))
        for (pos in 2:nrow(A)) {
          A[pos, pos - 1] = -1
        }
        beta_opt_k_list = vector(mode = "list", length = 2)
        W_opt_k = vector()
        # Monotonic increasing
        beta_opt_k = constrain_est(beta_hat = beta_hat_k, var_cov_hat = var_cov_hat_k, 
                                   A = A, increase = TRUE)
        # Discard taxa with structural zeros if testing the monotonic increasing order
        filter_ind = apply(zero_ind[[k]][, -1], 1, function(x) any(x == 1))
        beta_opt_k[filter_ind, ] = 0
        W_opt_k = cbind(W_opt_k, apply(beta_opt_k / sqrt(var_hat_k), 1, 
                                       function(x) abs(x)[length(x)]))
        beta_opt_k_list[[1]] = beta_opt_k
        # Monotonic decreasing
        beta_opt_k = constrain_est(beta_hat = beta_hat_k, var_cov_hat = var_cov_hat_k, 
                                   A = A, increase = FALSE)
        W_opt_k = cbind(W_opt_k, apply(beta_opt_k / sqrt(var_hat_k), 1, 
                                       function(x) abs(x)[length(x)]))
        beta_opt_k_list[[2]] = beta_opt_k
        # Which pattern does it present
        max_pos = apply(W_opt_k, 1, function(x) which.max(x))
        # Estimated coefficients
        beta_opt_k = vector()
        for (i in 1:n_taxa) {
          beta_opt_k = rbind(beta_opt_k, beta_opt_k_list[[max_pos[[i]]]][i, ])
        }
        beta_opt[[k]] = beta_opt_k
        # Test statistics under constraint
        W_opt_k = apply(W_opt_k, 1, function(x) max(x, na.rm = TRUE))
        W_opt[, k] = W_opt_k
        # Test statistics under null
        W_null_k = apply(beta_null_k / sqrt(var_hat_k[rep(seq(n_taxa), each = per_num), ]), 1, 
                         function(x) abs(x)[length(x)])
        W_null_k = matrix(W_null_k, nrow = n_taxa, ncol = per_num, byrow = TRUE)
        rownames(W_null_k) = taxa_id
        # P-values
        p_opt_k = apply(W_null_k - W_opt_k, 1, function(x) sum(x > 0) / per_num)
        names(p_opt_k) = taxa_id
        p_opt[, k] = p_opt_k
        # Q-values
        q_opt_k = p.adjust(p_opt_k, method = p_adj_method)
        names(q_opt_k) = taxa_id
        q_opt[, k] = q_opt_k
        # Nodes
        node_k = recode(max_pos, `1` = sum(k_ind), `2` = -sum(k_ind))
      } else if (pattern[k] == "tree") {
        # Tree order
        A = diag(1, nrow = sum(k_ind))
        beta_opt_k_list = vector(mode = "list", length = 2)
        W_opt_k = vector()
        # >= the reference group
        beta_opt_k = constrain_est(beta_hat = beta_hat_k, var_cov_hat = var_cov_hat_k, 
                                   A = A, increase = TRUE)
        # Discard taxa with structural zeros if testing >= the reference group
        filter_ind = apply(zero_ind[[k]][, -1], 1, function(x) any(x == 1))
        beta_opt_k[filter_ind, ] = 0
        W_opt_k = cbind(W_opt_k, apply(beta_opt_k / sqrt(var_hat_k), 1, 
                                       function(x) max(abs(x), na.rm = TRUE)))
        beta_opt_k_list[[1]] = beta_opt_k
        # <= the reference group
        beta_opt_k = constrain_est(beta_hat = beta_hat_k, var_cov_hat = var_cov_hat_k, 
                                   A = A, increase = FALSE)
        W_opt_k = cbind(W_opt_k, apply(beta_opt_k / sqrt(var_hat_k), 1, 
                                       function(x) max(abs(x), na.rm = TRUE)))
        beta_opt_k_list[[2]] = beta_opt_k
        # Which pattern does it present
        max_pos = apply(W_opt_k, 1, function(x) which.max(x))
        # Estimated coefficients
        beta_opt_k = vector()
        for (i in 1:n_taxa) {
          beta_opt_k = rbind(beta_opt_k, beta_opt_k_list[[max_pos[[i]]]][i, ])
        }
        beta_opt[[k]] = beta_opt_k
        # Test statistics under constraint
        W_opt_k = apply(W_opt_k, 1, function(x) max(x, na.rm = TRUE))
        W_opt[, k] = W_opt_k
        # Test statistics under null
        W_null_k = apply(beta_null_k / sqrt(var_hat_k[rep(seq(n_taxa), each = per_num), ]), 1, 
                         function(x) max(abs(x), na.rm = TRUE))
        W_null_k = matrix(W_null_k, nrow = n_taxa, ncol = per_num, byrow = TRUE)
        rownames(W_null_k) = taxa_id
        # P-values
        p_opt_k = apply(W_null_k - W_opt_k, 1, function(x) sum(x > 0) / per_num)
        names(p_opt_k) = taxa_id
        p_opt[, k] = p_opt_k
        # Q-values
        q_opt_k = p.adjust(p_opt_k, method = p_adj_method)
        names(q_opt_k) = taxa_id
        q_opt[, k] = q_opt_k
        # Nodes
        node_k = recode(max_pos, `1` = 1, `2` = -1)
      } else {
        # Umbrella order
        A0 = diag(1, nrow = sum(k_ind))
        for (pos in 2:nrow(A0)) {
          A0[pos, pos - 1] = -1
        }
        A = list()
        for (l in 1:sum(k_ind)) {
          A_l = A0
          A_l[-(1:l), ] = A_l[-(1:l), ] * (-1)
          A = c(A, list(A_l))
        }
        names(A) = paste0("A", 1:sum(k_ind))
        W_opt_k = vector()
        beta_opt_k_list = vector(mode = "list", length = sum(k_ind))
        # First increasing, then decreasing
        for (l in 1:sum(k_ind)) {
          # Estimated coefficients under constraint
          beta_opt_k = constrain_est(beta_hat = beta_hat_k, var_cov_hat = var_cov_hat_k, 
                                     A = A[[l]], increase = TRUE)
          # Discard taxa with structural zeros if the absolute abundance is increasing up to the node 
          filter_ind = apply(zero_ind[[k]][, 2:(l + 1), drop = FALSE], 1, function(x) any(x == 1))
          beta_opt_k[filter_ind, ] = 0
          beta_opt_k_list[[l]] = beta_opt_k
          pos1 = l; pos2 = sum(k_ind)
          W_opt_k = cbind(W_opt_k,
                          apply(beta_opt_k / sqrt(var_hat_k), 1, function(x) 
                            max(c(abs(x)[pos1], abs(x[pos1] - x[pos2])), na.rm = TRUE)))
        }
        # First decreasing, then increasing
        for (l in (sum(k_ind) + 1) : (2 * sum(k_ind))) {
          # Estimated coefficients under constraint
          beta_opt_k = constrain_est(beta_hat = beta_hat_k, var_cov_hat = var_cov_hat_k, 
                                     A = A[[l - sum(k_ind)]], increase = FALSE)
          beta_opt_k_list[[l]] = beta_opt_k
          pos1 = l - sum(k_ind); pos2 = sum(k_ind)
          W_opt_k = cbind(W_opt_k,
                          apply(beta_opt_k / sqrt(var_hat_k), 1, function(x) 
                            max(c(abs(x)[pos1], abs(x[pos1] - x[pos2])), na.rm = TRUE)))
        }
        # Which pattern does it present
        max_pos = apply(W_opt_k, 1, function(x) which.max(x))
        pos1 = sapply(max_pos, function(x) ifelse(x > sum(k_ind), x - sum(k_ind), x))
        pos2 = sum(k_ind)
        # Estimated coefficients
        beta_opt_k = vector()
        for (i in 1:n_taxa) {
          beta_opt_k = rbind(beta_opt_k, beta_opt_k_list[[max_pos[i]]][i, ])
        }
        beta_opt[[k]] = beta_opt_k
        # Test statistics under constraint
        W_opt_k = apply(W_opt_k, 1, function(x) max(x, na.rm = TRUE))
        W_opt[, k] = W_opt_k
        # Test statistics under null
        W_null_k = c()
        eff_size_k = beta_null_k / sqrt(var_hat_k[rep(seq(n_taxa), each = per_num), ])
        for (i in 1:n_taxa) {
          W_null_ik = eff_size_k[((i - 1) * per_num + 1):(i * per_num), ]
          W_null_ik = apply(W_null_ik, 1, function(x) 
            max(c(abs(x)[pos1[i]], abs(x[pos1[i]] - x[pos2])), na.rm = TRUE))
          W_null_k = c(W_null_k, W_null_ik)
        }
        W_null_k = matrix(W_null_k, nrow = n_taxa, ncol = per_num, byrow = TRUE)
        rownames(W_null_k) = taxa_id
        # P-values
        p_opt_k = apply(W_null_k - W_opt_k, 1, function(x) sum(x > 0) / per_num)
        names(p_opt_k) = taxa_id
        p_opt[, k] = p_opt_k
        # Q-values
        q_opt_k = p.adjust(p_opt_k, method = p_adj_method)
        names(q_opt_k) = taxa_id
        q_opt[, k] = q_opt_k
        # Nodes
        node_k = sapply(max_pos, function(x) ifelse(x > sum(k_ind), sum(k_ind) - x, x))
      }
      node_k[q_opt_k >= alpha] = 0
      node[, k] = node_k
    }
    beta_opt = Reduce('cbind', beta_opt)
    rownames(beta_opt) = taxa_id
    return(list(beta = beta_opt, W = W_opt, p_val = p_opt, q_val = q_opt, node = node))
  }
  
  # Trend test summary
  if (!is.null(pattern)) {
    if (is.null(group)) stop("Trend test requires specifying the group variables.") 
    res_pattern = trend_test(para = group, pattern = pattern)  
    res_pattern$diff_abn = ifelse(res_pattern$q_val < alpha, TRUE, FALSE)
  } else {
    res_pattern = NULL
  }
  
  # 5. Combine the information of structural zeros
  # Set p/q-values of structural zeros to be 0s.
  if (struc_zero) {
    for (k in 1:length(group_all)) {
      k_ind = grepl(group_all[k], covariates)
      zero_mask_k = 1 - apply(zero_ind[[k]], 1, function(x) any(x == 1))
      res$p_val[, k_ind] = res$p_val[, k_ind] * zero_mask_k
      res$q_val[, k_ind] = res$q_val[, k_ind] * zero_mask_k
    }
    res$diff_abn = ifelse(res$q_val < alpha, TRUE, FALSE)
    
    # Global test
    if (global) {
      for (k in 1:length(group)) {
        zero_mask_k = 1 - apply(zero_ind[[k]], 1, function(x) any(x == 1))
        res_global$p_val[, k] = res_global$p_val[, k] * zero_mask_k
        res_global$q_val[, k] = res_global$q_val[, k] * zero_mask_k
      }
      res_global$diff_abn = ifelse(res_global$q_val < alpha, TRUE, FALSE)
    }
    
    # Directional test
    if (direct) {
      if (dunnett) {
        for (k in 1:length(group)) {
          k_ind = grepl(group[k], colnames(res_direct$beta))
          zero_mask_k = 1 - (zero_ind[[k]][, -1] - zero_ind[[k]][, 1]) # equals to 0, 1, 2
          zero_mask_k[zero_mask_k == 2] = 0 
          res_direct$p_val[, k_ind] = res_direct$p_val[, k_ind] * zero_mask_k
          res_direct$q_val[, k_ind] = res_direct$q_val[, k_ind] * zero_mask_k
        }
        res_direct$diff_abn = ifelse(res_direct$q_val < alpha, TRUE, FALSE)
      } else {
        for (k in 1:length(group)) {
          k_ind = grepl(group[k], colnames(res_direct$beta))
          zero_mask_k = t(apply(zero_ind[[k]], 1, function(x) combn_fun(x, fun = diff)))
          zero_mask_k = zero_mask_k[, -(1:ncol(zero_ind[[k]]))]
          zero_mask_k = 1 - zero_mask_k
          zero_mask_k[zero_mask_k ==2] = 0
          res_direct$p_val[, k_ind] = res_direct$p_val[, k_ind] * zero_mask_k
          res_direct$q_val[, k_ind] = res_direct$q_val[, k_ind] * zero_mask_k
        }
        res_direct$diff_abn = ifelse(res_direct$q_val < alpha, TRUE, FALSE)
      }
    }
  }

  # 6. Outputs
  out = list(feature_table = feature_table, zero_ind = zero_ind_all, samp_frac = d_hat, 
             delta_em = delta_em, delta_wls = delta_wls, res = res, 
             res_global = res_global, res_direct = res_direct, res_pattern = res_pattern)
  return(out)
}


