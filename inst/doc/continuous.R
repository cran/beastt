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
set.seed(1234)

summary(int_norm_df)
summary(ex_norm_df)
sd_external_control <- 0.15
sd_internal_control <- 0.15
sd_internal_treated <- 0.15

## -----------------------------------------------------------------------------
ps_model <- ~ cov1 + cov2 + cov3 + cov4
ps_obj <- calc_prop_scr(internal_df = filter(int_norm_df, trt == 0), 
                        external_df = ex_norm_df, 
                        id_col = subjid,
                        model = ps_model)

ps_obj

## -----------------------------------------------------------------------------
prop_scr_hist(ps_obj)
prop_scr_dens(ps_obj, variable = "ipw")

## -----------------------------------------------------------------------------
prop_scr_love(ps_obj, reference_line = 0.1)

## -----------------------------------------------------------------------------
pwr_prior <- calc_power_prior_norm(ps_obj,
                                   response = y,
                                   prior = dist_normal(0.5, 10), 
                                   external_sd = sd_external_control)
plot_dist(pwr_prior)

## -----------------------------------------------------------------------------
n_external <- nrow(ex_norm_df)
mix_prior <- robustify_norm(pwr_prior, n_external, weights = c(0.5, 0.5))
plot_dist("Power Prior" = pwr_prior,
          "Vague Prior" = dist_normal(mean = mix_means(mix_prior)["vague"],
                                      sd = mix_sigmas(mix_prior)["vague"]),
          "Robust Mixture Prior" = mix_prior)

## -----------------------------------------------------------------------------
post_control <- calc_post_norm(filter(int_norm_df, trt == 0),
                               response = y, 
                               prior = mix_prior,
                               internal_sd = sd_internal_control)
plot_dist(post_control)

## -----------------------------------------------------------------------------
post_treated <- calc_post_norm(internal_data = filter(int_norm_df, trt == 1),
                               response = y,
                               prior = dist_normal(mean = mix_means(mix_prior)["vague"],
                                                   sd = mix_sigmas(mix_prior)["vague"]),
                               internal_sd = sd_internal_treated)
plot_dist("Control Posterior" = post_control,
          "Treatment Posterior" = post_treated)

## -----------------------------------------------------------------------------
c(mean = mean(post_control),
  median = median(post_control),
  variance = variance(post_control))

## -----------------------------------------------------------------------------
hdr(post_control)        # 95% HDR
hdr(post_control, 90)    # 90% HDR

## -----------------------------------------------------------------------------
hilo(post_control)       # 95% credible interval
hilo(post_control, 90)   # 90% credible interval
quantile(post_control, c(.025, .975))[[1]]    # 95% CrI via quantile function

## -----------------------------------------------------------------------------
cdf(post_control, q = .7)   # Pr(theta_C < 0.7 | D)

## -----------------------------------------------------------------------------
density(post_control, at = .7)                # density at 0.7
density(post_control, at = .7, log = TRUE)    # log density at 0.7

## -----------------------------------------------------------------------------
samp_control <- generate(x = post_control, times = 100000)[[1]]
ggplot(data.frame(samp = samp_control), aes(x = samp)) +
  labs(y = "Density", x = expression(theta[C])) +
  ggtitle(expression(paste("Posterior Samples of ", theta[C]))) +
  geom_histogram(aes(y = after_stat(density)), color = "#5398BE", fill = "#5398BE",
                 position = "identity", binwidth = .01, alpha = 0.5) +
  geom_density(color = "black") +
  theme_bw()

## -----------------------------------------------------------------------------
samp_treated <- generate(x = post_treated, times = 100000)[[1]]
ggplot(data.frame(samp = samp_treated), aes(x = samp)) +
  labs(y = "Density", x = expression(theta[T])) +
  ggtitle(expression(paste("Posterior Samples of ", theta[T]))) +
  geom_histogram(aes(y = after_stat(density)), color = "#FFA21F", fill = "#FFA21F",
                 position = "identity", binwidth = .01, alpha = 0.5) +
  geom_density(color = "black") +
  theme_bw()

## -----------------------------------------------------------------------------
samp_trt_diff <- samp_treated - samp_control
ggplot(data.frame(samp = samp_trt_diff), aes(x = samp)) +
  labs(y = "Density", x = expression(paste(theta[T], " - ", theta[C]))) +
  ggtitle(expression(paste("Posterior Samples of ", theta[T], " - ", theta[C]))) +
  geom_histogram(aes(y = after_stat(density)), color = "#FF0000", fill = "#FF0000",
                 position = "identity", binwidth = .01, alpha = 0.5) +
  geom_density(color = "black") +
  theme_bw()

## -----------------------------------------------------------------------------
mean(samp_trt_diff > 0)

## -----------------------------------------------------------------------------
parameters(post_treated)

## -----------------------------------------------------------------------------
mix_means(post_control)    # means of each normal component
mix_sigmas(post_control)   # SDs of each normal component

## -----------------------------------------------------------------------------
parameters(post_control)$w[[1]]

## -----------------------------------------------------------------------------
post_control_no_brrw <- calc_post_norm(filter(int_norm_df, trt == 0),
                                       response = y,
                                       prior = dist_normal(mean = mix_means(mix_prior)["vague"],
                                                           sd = mix_sigmas(mix_prior)["vague"]),
                                       internal_sd = sd_internal_control)
n_int_ctrl <- nrow(filter(int_norm_df, trt == 0))   # sample size of internal control arm
var_no_brrw <- variance(post_control_no_brrw)       # post variance of theta_C without borrowing
var_brrw <- variance(post_control)                  # post variance of theta_C with borrowing
ess <- n_int_ctrl * var_no_brrw / var_brrw          # effective sample size
ess

