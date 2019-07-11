rm(list=ls())
# code to fit models for fragmentation synthesis
# so far: bayesian framework for approximately ML-like results (i.e.,
# with non-informative priors - against the science of Gelman)
library(tidyverse)
library(brms)

# load the data
frag <- read_csv(paste0(path2data, '1_biodiv_frag_fcont_10_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))


#----- simplest model: diversity as a function of fragment size; allow fragment size to vary by study----
get_prior(S_std_2 ~ c.lfs + (c.lfs | dataset_label),
          family = lognormal(),
          data = frag)

# set some weakly regularising priors
# rp <- c(prior(normal(0, 1), class = Intercept),
#         prior(normal(0, 0.5), class = b),
#         prior(exponential(1), class = sd))
# 

Sstd2_lognorm_fragSize <- brm(S_std_2 ~ c.lfs + (c.lfs | dataset_label), 
                              # fit to data with variation in frag_size_num
                              data = frag,
                              #prior = rp,
                              family = 'lognormal', # our standardised richness are not integer values
                              cores = 4, chains = 4)

# refit without the studies with pooled sampling designs
Sstd2_lognorm_fragSize_pool <- brm(S_std_2 ~ c.lfs + (c.lfs | dataset_label), 
                              # fit to data with variation in frag_size_num
                              data = frag %>% filter(sample_design!='pooled'),
                              #prior = rp,
                              family = 'lognormal', # our standardised richness are not integer values
                              cores = 4, chains = 4)

Sn_lognorm_fragSize <- brm(S_n ~ c.lfs + (c.lfs | dataset_label), 
                           data = frag,
                           #prior = rp,
                           family = 'lognormal',
                           cores = 4, chains = 4)
# refit without the studies with pooled sampling designs
Sn_lognorm_fragSize_pool <- brm(S_n ~ c.lfs + (c.lfs | dataset_label), 
                           data = frag %>% filter(sample_design!='pooled'),
                           #prior = rp,
                           family = 'lognormal',
                           cores = 4, chains = 4)

Scov_lognorm_fragSize <- brm(S_cov ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag,
                             #prior = rp,
                             family = hurdle_lognormal(),
                             cores = 4, chains = 4)

# refit without studies with pooled sampling designs
Scov_lognorm_fragSize_pool <- brm(S_cov ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag %>% filter(sample_design!='pooled'),
                             #prior = rp,
                             family = hurdle_lognormal(),
                             cores = 4, chains = 4)


S_PIE_lognorm_fragSize <- brm(S_PIE ~ c.lfs + (c.lfs | dataset_label), 
                              data = frag,
                              #prior = rp,
                              family = 'lognormal',
                              cores = 4, chains = 4)

# refit without studies with pooled sampling designs
S_PIE_lognorm_fragSize_pool <- brm(S_PIE ~ c.lfs + (c.lfs | dataset_label), 
                              data = frag %>% filter(sample_design!='pooled'),
                              #prior = rp,
                              family = 'lognormal',
                              cores = 4, chains = 4)

S_chao_lognorm_fragSize <- brm(S_chao ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag,
                               #prior = rp,
                               family = 'lognormal',
                               cores = 4, chains = 4)

# refit without the studies with pooled sampling designs
S_chao_lognorm_fragSize_pool <- brm(S_chao ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag %>% filter(sample_design!='pooled'),
                               #prior = rp,
                               family = 'lognormal',
                               cores = 4, chains = 4)


Nstd_lognorm_fragSize <- brm(N_std ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag,
                             #prior = rp,
                             family = 'lognormal',
                             cores = 4, chains = 4)


Nstd_lognorm_fragSize_pool <- brm(N_std ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag %>% filter(sample_design!='pooled'),
                             #prior = rp,
                             family = 'lognormal',
                             cores = 4, chains = 4)


save(Sstd2_lognorm_fragSize,
     Sstd2_lognorm_fragSize_pool,
     Sn_lognorm_fragSize,
     Sn_lognorm_fragSize_pool,
     S_PIE_lognorm_fragSize,
     S_PIE_lognorm_fragSize_pool,
     Scov_lognorm_fragSize,
     Scov_lognorm_fragSize_pool,
     S_chao_lognorm_fragSize,
     S_chao_lognorm_fragSize_pool,
     Nstd_lognorm_fragSize,
     Nstd_lognorm_fragSize_pool,
     file = '~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms_ref_new.Rdata')
