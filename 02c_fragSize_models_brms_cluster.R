rm(list=ls())
# code to fit models for fragmentation synthesis
# so far: bayesian framework for approximately ML-like results (i.e.,
# with non-informative priors - against the science of Gelman)
library(tidyverse)
library(brms)

# load the data
frag <- read_csv('~/Dropbox/Frag Database (new)/files_datapaper/Analysis/2_biodiv_frag_fcont_10_mabund_as_is.csv')
# cluster version
# frag <- read_csv('/gpfs1/data/idiv_chase/sablowes/fragmentation/data/2_biodiv_frag_fcont_10_mabund_as_is.csv')

##--create some covariates for easier workflow--
# mean-centred log(fragment.size)
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))
# standardised log(fragment.size) may improve sampling if problems arise...
# frag$slfs <- scale(log(frag$frag_size_num),center = T, scale = T)
# add log-transformed (to fit multivariate normal?)
# frag$lSstd1 <- log(frag$S_std_1)
# frag$lSstd2 <- log(frag$S_std_2)
# frag$lSn <- log(frag$S_n)
# frag$lScov <- log(frag$S_cov)
# frag$lS_PIE <- log(frag$S_PIE)
# frag$lS_chao <- log(frag$S_chao)
# frag$lNstd <- log(frag$N_std)


#----- simplest model: diversity as a function of fragment size; allow fragment size to vary by study----
get_prior(S_std_2 ~ c.lfs + (c.lfs | dataset_label),
          family = lognormal(),
          data = frag)

with(frag, lm(log(Sn) ~ c.lfs))

# set some weakly regularising priors
rp <- c(prior(normal(2,2), class = Intercept),
        prior(normal(0,1), class = b),
        prior(exponential(1), class = sd))

# not sure which of the standardised richness values to model: Felix say Sstd2
# Sstd1_lognorm_fragSize <- brm(S_std_1 ~ c.lfs + (c.lfs | dataset_label), 
#                               # fit to data with variation in frag_size_num
#                               data = frag,
#                               prior = rp,
#                               family = 'lognormal',
#                               cores = 4, chains = 4)
# 
# pp_check(Sstd1_lognorm_fragSize) +
#   scale_x_continuous(trans = 'log', 
#                      breaks = c(1,8,16,32,64,128))

Sstd2_lognorm_fragSize <- brm(S_std_2 ~ c.lfs + (c.lfs | dataset_label), 
                              # fit to data with variation in frag_size_num
                              data = frag,
                              prior = rp,
                              family = 'lognormal', # our standardised richness are not integer values
                              cores = 4, chains = 4)

plot(Sstd2_lognorm_fragSize) # intercept(s) don't sample great, but we look okish based on eff. sample sizes
                             # there was a small improvement with prior_intercept ~ norm(2, 2) cf. norm(0,2)
pp_check(Sstd2_lognorm_fragSize) +
  scale_x_continuous(trans = 'log')

Sn_lognorm_fragSize <- brm(S_n ~ c.lfs + (c.lfs | dataset_label), 
                           data = frag,
                           prior = rp,
                           family = 'lognormal',
                           cores = 4, chains = 4)

plot(Sn_lognorm_fragSize) # intercept(s) don't sample great, but we look okish based on eff. sample sizes
pp_check(Sn_lognorm_fragSize) +
  scale_x_continuous(trans = 'log') # smooth :)

Scov_lognorm_fragSize <- brm(S_cov ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag,
                             prior = rp,
                             family = hurdle_lognormal(),
                             cores = 4, chains = 4)

plot(Scov_lognorm_fragSize) # intercept(s) don't sample great, but we look okish based on eff. sample sizes
pp_check(Sn_lognorm_fragSize) +
  scale_x_continuous(trans = 'log') # smooth :)

S_PIE_lognorm_fragSize <- brm(S_PIE ~ c.lfs + (c.lfs | dataset_label), 
                              data = frag,
                              prior = rp,
                              family = 'lognormal',
                              cores = 4, chains = 4)

plot(S_PIE_lognorm_fragSize) # intercept(s) don't sample great, but we look okish based on eff. sample sizes
pp_check(S_PIE_lognorm_fragSize) +
  scale_x_continuous(trans = 'log') # smooth :)

S_chao_lognorm_fragSize <- brm(S_chao ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag,
                               prior = rp,
                               family = 'lognormal',
                               cores = 4, chains = 4)

plot(S_chao_lognorm_fragSize) # intercept least good, but eff samples > 200
pp_check(S_chao_lognorm_fragSize) +
  scale_x_continuous(trans = 'log') # good

Nstd_lognorm_fragSize <- brm(N_std ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag,
                             prior = rp,
                             family = 'lognormal',
                             cores = 4, chains = 4)

plot(Nstd_lognorm_fragSize) # increasing the mean of the prior for the intercept further could be warranted?
pp_check(Nstd_lognorm_fragSize) +
  scale_x_continuous(trans = 'log') # but I think we are ok to proceed

save(#Sstd1_lognorm_fragSize, 
     Sstd2_lognorm_fragSize,
     Sn_lognorm_fragSize,
     S_PIE_lognorm_fragSize,
     Scov_lognorm_fragSize,
     S_chao_lognorm_fragSize,
     Nstd_lognorm_fragSize,
     file = '~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms_ref.Rdata')
