library(tidyverse)

n_taxa = 10
n_sample = 100

set.seed(1)
alpha_0 = sample(c(log(50), log(200), log(10000)), 
                 n_taxa, replace = TRUE, 
                 prob = c(0.6, 0.3, 0.1)) 

alpha_1 = rnorm(n_taxa, mean = 0, sd = 3)
alpha_2 = rnorm(n_taxa, mean = 0, sd = 2)

meta_data = data.frame(x1 = rep(c(0, 1), each = n_sample/2),
                       x2 = runif(n_sample, 10, 20))

# Case 1: only x1
A = matrix(NA, nrow = n_taxa, ncol = n_sample)

for (i in 1:n_taxa) {
  A[i, ] = alpha_0[i] + alpha_1[i] * meta_data$x1 
}

A_alr = t(t(A) - A[n_taxa, ])
A_alr = A_alr[-n_taxa, ]

beta = apply(A_alr, 1, function(x) coef(lm(x ~ meta_data$x1))[2])

dat_fig = data.frame(sample_id = factor(1:9), 
                     alpha_1 = rank(alpha_1[-n_taxa]), beta_1 = rank(beta))

p1 = dat_fig %>% ggplot(aes(x = alpha_1, y = beta_1)) + 
  geom_point(aes(color = sample_id)) +
  labs(x = "Alpha1", y = "Beta1") +
  scale_color_discrete(name = "Sample") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw()

p1

# Case 2: both x1 and x2
A = matrix(NA, nrow = n_taxa, ncol = n_sample)

for (i in 1:n_taxa) {
  A[i, ] = alpha_0[i] + alpha_1[i] * meta_data$x1 + alpha_2[i] * meta_data$x2
}

A_alr = t(t(A) - A[n_taxa, ])
A_alr = A_alr[-n_taxa, ]

beta = apply(A_alr, 1, function(x) coef(lm(x ~ meta_data$x1 + meta_data$x2))[2])

dat_fig = data.frame(sample_id = factor(1:9), 
                     alpha_1 = rank(alpha_1[-n_taxa]), beta_1 = rank(beta))

p2 = dat_fig %>% ggplot(aes(x = alpha_1, y = beta_1)) + 
  geom_point(aes(color = sample_id)) +
  labs(x = "Alpha1", y = "Beta1") +
  scale_color_discrete(name = "Sample") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw()

p2

