## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(beastt)

## ----class.source = 'fold-hide'-----------------------------------------------
library(tibble)
library(distributional)
library(dplyr)
library(ggplot2)
library(rstan)
set.seed(1234)

summary(int_tte_df)
summary(ex_tte_df)

## -----------------------------------------------------------------------------
ps_obj <- calc_prop_scr(internal_df = filter(int_tte_df, trt == 0),
                        external_df = ex_tte_df,
                        id_col = subjid,
                        model = ~ cov1 + cov2 + cov3 + cov4)
ps_obj

## -----------------------------------------------------------------------------
prop_scr_hist(ps_obj)
prop_scr_dens(ps_obj, variable = "ipw")

## -----------------------------------------------------------------------------
prop_scr_love(ps_obj, reference_line = 0.1)

## -----------------------------------------------------------------------------
pwr_prior <- calc_power_prior_weibull(ps_obj,
                                      response = y,
                                      event = event,
                                      intercept = dist_normal(0, 10),
                                      shape = 50,
                                      approximation = "Laplace")
plot_dist(pwr_prior)

## -----------------------------------------------------------------------------
r_external <- sum(ex_tte_df$event)   # number of observed events
mix_prior <- robustify_mvnorm(pwr_prior, r_external, weights = c(0.5, 0.5))   # RMP
mix_means(mix_prior)     # mean vectors
mix_sigmas(mix_prior)    # mean covariance matrices
#plot_dist(mix_prior)

## -----------------------------------------------------------------------------
post_control <- calc_post_weibull(filter(int_tte_df, trt == 0),
                                  response = y,
                                  event = event,
                                  prior = mix_prior,
                                  analysis_time = 12)
summary(post_control)$summary
#plot_dist(post_control)

## -----------------------------------------------------------------------------
surv_prob_control <- as.data.frame(extract(post_control, pars = c("survProb")))[,1]
ggplot(data.frame(samp = surv_prob_control), aes(x = samp)) +
  labs(y = "Density", x = expression(paste(S[C], "(t=12)"))) +
  ggtitle(expression(paste("Posterior Samples of ", S[C], "(t=12)"))) +
  geom_histogram(aes(y = after_stat(density)), color = "#5398BE", fill = "#5398BE",
                 position = "identity", binwidth = .01, alpha = 0.5) +
  geom_density(color = "black") +
  coord_cartesian(xlim = c(-0.2, 0.8)) +
  theme_bw()

## -----------------------------------------------------------------------------
vague_prior <- dist_multivariate_normal(mu = list(mix_means(mix_prior)[[2]]),
                                        sigma = list(mix_sigmas(mix_prior)[[2]]))
post_treated <- calc_post_weibull(filter(int_tte_df, trt == 1),
                                  response = y,
                                  event = event,
                                  prior = vague_prior,
                                  analysis_time = 12)
summary(post_treated)$summary
#plot_dist(post_treated)

## -----------------------------------------------------------------------------
surv_prob_treated <- as.data.frame(extract(post_treated, pars = c("survProb")))[,1]
ggplot(data.frame(samp = surv_prob_treated), aes(x = samp)) +
  labs(y = "Density", x = expression(paste(S[T], "(t=12)"))) +
  ggtitle(expression(paste("Posterior Samples of ", S[T], "(t=12)"))) +
  geom_histogram(aes(y = after_stat(density)), color = "#FFA21F", fill = "#FFA21F",
                 position = "identity", binwidth = .01, alpha = 0.5) +
  geom_density(color = "black") +
  coord_cartesian(xlim = c(-0.2, 0.8)) +
  theme_bw()

## -----------------------------------------------------------------------------
samp_trt_diff <- surv_prob_treated - surv_prob_control
ggplot(data.frame(samp = samp_trt_diff), aes(x = samp)) +
  labs(y = "Density", x = expression(paste(S[T], "(t=12) - ", S[C], "(t=12)"))) +
  ggtitle(expression(paste("Posterior Samples of ", S[T],
                           "(t=12) - ", S[C], "(t=12)"))) +
  geom_histogram(aes(y = after_stat(density)), color = "#FF0000", fill = "#FF0000",
                 position = "identity", binwidth = .01, alpha = 0.5) +
  geom_density(color = "black") +
  coord_cartesian(xlim = c(-0.2, 0.8)) +
  theme_bw()

## -----------------------------------------------------------------------------
mean(samp_trt_diff > 0)

## -----------------------------------------------------------------------------
c(mean = mean(samp_trt_diff),
  median = median(samp_trt_diff),
  SD = sd(samp_trt_diff))

## -----------------------------------------------------------------------------
quantile(samp_trt_diff, c(.025, .975))    # 95% CrI

## -----------------------------------------------------------------------------
post_ctrl_no_brrw <- calc_post_weibull(filter(int_tte_df, trt == 0),
                                       response = y,
                                       event = event,
                                       prior = vague_prior,
                                       analysis_time = 12)
surv_prob_ctrl_nb <- as.data.frame(extract(post_ctrl_no_brrw, pars = c("survProb")))[,1]
n_int_ctrl <- nrow(filter(int_tte_df, trt == 0))  # sample size of internal control arm
var_no_brrw <- var(surv_prob_ctrl_nb)             # post variance of S_C(t) without borrowing
var_brrw <- var(surv_prob_control)                # post variance of S_C(t) with borrowing
ess <- n_int_ctrl * var_no_brrw / var_brrw        # effective sample size
ess

