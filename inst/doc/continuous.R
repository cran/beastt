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
n_external <- nrow(ex_norm_df)

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
                                   prior = dist_normal(50, 10), 
                                   external_sd = sd_external_control)
plot_dist(pwr_prior)

## -----------------------------------------------------------------------------
post <- calc_post_norm(ps_obj,
                       response = y, 
                       prior = pwr_prior,
                       internal_sd = sd_internal_control)
plot_dist(post)

## -----------------------------------------------------------------------------
mixed_prior <- robustify_norm(pwr_prior, n_external)

post_mixed <- calc_post_norm(ps_obj,
                             response = y, 
                             prior= mixed_prior,
                             internal_sd = sd_internal_control)
plot_dist("Control Posterior" = post, 
          "Mixed Posterior" = post_mixed)

## -----------------------------------------------------------------------------
sd_internal_treated <- 0.15
post_treated <- calc_post_norm(internal_data = filter(int_norm_df, trt == 1),
                               response = y,
                               prior = dist_normal(50, 10),
                               internal_sd = sd_internal_treated)

plot_dist("Control Posterior" = post, 
          "Mixed Posterior" = post_mixed,
          "Treatment Posterior" = post_treated)

