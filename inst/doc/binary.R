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
post <- calc_post_beta(ps_obj,
                       response = y,
                       prior = pwr_prior)
plot_dist(post)

## -----------------------------------------------------------------------------
mix_prior <- dist_mixture(pwr_prior,
                          dist_beta(0.5, 0.5),
                          weights = c(0.5, 0.5))
post_mixed <- calc_post_beta(ps_obj,
                             response = y,
                             prior = mix_prior)
plot_dist("Control Posterior" = post, 
          "Mixed Posterior" = post_mixed)

## -----------------------------------------------------------------------------
post_treated <- calc_post_beta(internal_data = filter(int_binary_df, trt == 1),
                               response = y,
                               prior = dist_beta(0.5, 0.5))
plot_dist("Control Posterior" = post, 
          "Mixed Posterior" = post_mixed,
          "Treatment Posterior" = post_treated)

