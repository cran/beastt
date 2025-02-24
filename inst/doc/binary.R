## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(beastt)
library(distributional)
library(dplyr)
library(ggplot2)
set.seed(1234)

summary(int_binary_df)
summary(ex_binary_df)

## -----------------------------------------------------------------------------
ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                        external_df = ex_binary_df,
                        id_col = subjid,
                        model = ~ cov1 + cov2 + cov3 + cov4)
ps_obj

## -----------------------------------------------------------------------------
prop_scr_hist(ps_obj)
prop_scr_dens(ps_obj, variable = "ipw")

## -----------------------------------------------------------------------------
prop_scr_love(ps_obj, reference_line = 0.1)

## -----------------------------------------------------------------------------
pwr_prior <- calc_power_prior_beta(ps_obj,
                                   response = y,
                                   prior = dist_beta(0.5, 0.5))
plot_dist(pwr_prior)

## -----------------------------------------------------------------------------
vague_prior <- dist_beta(0.5, 0.5)
mix_prior <- dist_mixture(informative = pwr_prior,
                          vague = vague_prior,
                          weights = c(0.5, 0.5))
plot_dist("Power Prior" = pwr_prior,
          "Vague Prior" = vague_prior,
          "Robust Mixture Prior" = mix_prior)

## -----------------------------------------------------------------------------
post_control <- calc_post_beta(filter(int_binary_df, trt == 0),
                               response = y,
                               prior = mix_prior)
plot_dist(post_control)

## -----------------------------------------------------------------------------
post_treated <- calc_post_beta(internal_data = filter(int_binary_df, trt == 1),
                               response = y,
                               prior = vague_prior)
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
cdf(post_control, q = .5)   # Pr(theta_C < 0.5 | D)

## -----------------------------------------------------------------------------
density(post_control, at = .5)                # density at 0.5
density(post_control, at = .5, log = TRUE)    # log density at 0.5

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
parameters(post_control)$w[[1]]

## -----------------------------------------------------------------------------
post_control_no_brrw <- calc_post_beta(filter(int_binary_df, trt == 0),
                                       response = y,
                                       prior = vague_prior)
n_int_ctrl <- nrow(filter(int_binary_df, trt == 0))   # sample size of internal control arm
var_no_brrw <- variance(post_control_no_brrw)    # post variance of theta_C without borrowing
var_brrw <- variance(post_control)               # post variance of theta_C with borrowing
ess <- n_int_ctrl * var_no_brrw / var_brrw       # effective sample size
ess

