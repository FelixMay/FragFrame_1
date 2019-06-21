rm(list=ls())
# code to fit models for fragmentation synthesis
# so far: bayesian framework for approximately ML-like results (i.e.,
# with non-informative priors - against the science of Gelman)
library(tidyverse)
library(brms)

# load the data
frag <- read_csv('~/Dropbox/Frag Database (new)/files_datapaper/Analysis/2_biodiv_frag_fcont_10_mabund_as_is.csv')

# load the meta data
meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')

##--create some covariates for easier workflow--
# mean-centred log(fragment.size)
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))
# standardised log(fragment.size) may improve sampling if problems arise...
frag$slfs <- scale(log(frag$frag_size_num),center = T, scale = T)
# add log-transformed (to fit multivariate normal?)
frag$lSstd1 <- log(frag$S_std_1)
frag$lSstd2 <- log(frag$S_std_2)
frag$lSn <- log(frag$S_n)
frag$lScov <- log(frag$S_cov)
frag$lS_PIE <- log(frag$S_PIE)
frag$lS_chao <- log(frag$S_chao)
frag$lNstd <- log(frag$N_std)

# how steep is the slope? Think about regularising prior (plot for completeness)
with(frag, summary(lm(lSstd2~c.lfs)))
with(frag, summary(lm(lSstd2~c.lfs*continent7)))
with(frag, car::Anova(lm(lSstd2~c.lfs*continent7)))
with(frag, car::Anova(lm(lSstd2~c.lfs*continent8)))
with(frag, car::Anova(lm(lSstd2~c.lfs*climate)))
with(frag, AIC(lm(lSstd1~c.lfs*continent7)))
with(frag, AIC(lm(lSstd1~c.lfs*continent8))) # deltaAIC ~ -7: support for central and north america being separate
with(frag, AIC(lm(lSstd1~c.lfs*climate))) # oversimplification: deltaAIC ~ +40

# how many studies per continent8
frag %>% 
  group_by(continent8) %>% 
  summarise(n_study = n_distinct(dataset_label))

# create new grouping covariates: concatenate continent8_matrix_time_taxa_biome
frag <- frag %>% 
  unite(cbmtt, c(continent8, biome, Matrix.category, time.since.fragmentation, taxa),
        remove = F) %>% 
  unite(bmtt, c(biome, Matrix.category, time.since.fragmentation, taxa),
        remove = F)

#----- simplest model: diversity as a function of fragment size; allow fragment size to vary by study----
get_prior(S_std_2 ~ c.lfs + (c.lfs | dataset_label),
          family = lognormal(),
          data = frag)
# set some weakly regularising priors
rp <- c(prior(normal(0,2), class = Intercept),
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
                              family = 'lognormal',
                              cores = 4, chains = 4)

load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms.Rdata')
load('~/Dropbox/1current/fragmentation_synthesis/results/Sstd2_model_explore.Rdata')

Sstd2_lognorm_fS_climate <- update(Sstd2_lognorm_fragSize,
                                 formula. = ~ c.lfs*climate + (c.lfs*climate | dataset_label), 
                                 newdata = frag,
                                 # prior = rp,
                                 family = 'lognormal',
                                 cores = 4, chains = 4)

Sstd2_lognorm_fS_clim_cont<- update(Sstd2_lognorm_fragSize,
                                 formula. = ~ c.lfs*climate + (c.lfs*climate | continent8/bmtt/dataset_label), 
                                 newdata = frag,
                                 # prior = rp,
                                 family = 'lognormal',
                                 cores = 4, chains = 4)

Sstd2_lognorm_fS_clim_cont_bmtt <- update(Sstd2_lognorm_fragSize,
                                  formula. = ~ c.lfs + (c.lfs | climate/continent8/bmtt/dataset_label), 
                                  newdata = frag,
                                  # prior = rp,
                                  family = 'lognormal',
                                  cores = 4, chains = 4)

Sstd2_lognorm_fS_cbmtt0 <- update(Sstd2_lognorm_fragSize,
                                     formula. = ~ 1 + (c.lfs | cbmtt/dataset_label), 
                                     newdata = frag,
                                     # prior = rp,
                                     family = 'lognormal',
                                     cores = 4, chains = 4)

# this one makes predictions comparable to others, but does not draw many effective samples
Sstd2_lognorm_fS_continent0 <- update(Sstd2_lognorm_fragSize,
                                      formula. = ~ 1 + (c.lfs | continent8/dataset_label), 
                                      newdata = frag,
                                      # prior = rp,
                                      family = 'lognormal',
                                      cores = 4, chains = 4)


Sstd2_lognorm_fS_continent1 <- update(Sstd2_lognorm_fragSize,
                                    formula. = ~ c.lfs + (c.lfs | continent8/dataset_label), 
                                    newdata = frag,
                                    # prior = rp,
                                    family = 'lognormal',
                                    cores = 4, chains = 4)

Sstd2_lognorm_fS_continent_bmtt <- update(Sstd2_lognorm_fragSize,
                                     formula. = ~ c.lfs + (c.lfs | continent8/bmtt/dataset_label), 
                                     newdata = frag,
                                     # prior = rp,
                                     family = 'lognormal',
                                     cores = 4, chains = 4)

Sstd2_lognorm_fS_cbmtt1 <- update(Sstd2_lognorm_fragSize,
                                      formula. = ~ c.lfs + (c.lfs | cbmtt/dataset_label), 
                                      newdata = frag,
                                      # prior = rp,
                                      family = 'lognormal',
                                      cores = 4, chains = 4)

Sstd2_lognorm_fS_clim_cbmtt0 <- update(Sstd2_lognorm_fragSize,
                                       formula. = ~ 1 + (c.lfs | climate/cbmtt/dataset_label), 
                                       newdata = frag,
                                       # prior = rp,
                                       family = 'lognormal',
                                       cores = 4, chains = 4)
# refit with default priors to see if sampling improves
Sstd2_lognorm_fS_clim_cbmtt0 <- brm(S_std_2 ~ 1 + (c.lfs | climate/cbmtt/dataset_label), 
                                       data = frag,
                                       # prior = rp,
                                       family = 'lognormal',
                                       cores = 4, chains = 4)

Sstd2_lognorm_fS_clim_cbmtt1 <- update(Sstd2_lognorm_fragSize,
                                      formula. = ~ c.lfs + (c.lfs | climate/cbmtt/dataset_label), 
                                      newdata = frag,
                                      # prior = rp,
                                      family = 'lognormal',
                                      cores = 4, chains = 4)
brms::model_weights(Sstd2_lognorm_fragSize,
     # Sstd2_lognorm_fS_clim_cont,
     # Sstd2_lognorm_fS_clim_cont_bmtt,
     # Sstd2_lognorm_fS_cbmtt0,
     # Sstd2_lognorm_fS_continent0,
     Sstd2_lognorm_fS_continent1,
     Sstd2_lognorm_fS_continent_bmtt,
     Sstd2_lognorm_fS_cbmtt1,
     # Sstd2_lognorm_fS_clim_cbmtt0,
     # Sstd2_lognorm_fS_clim_cbmtt1, 
     weights = 'waic') 

brms::waic(Sstd2_lognorm_fragSize,
                    # Sstd2_lognorm_fS_clim_cont,
                    # Sstd2_lognorm_fS_clim_cont_bmtt,
                    # Sstd2_lognorm_fS_cbmtt0,
                    # Sstd2_lognorm_fS_continent0,
                    Sstd2_lognorm_fS_continent1,
                    Sstd2_lognorm_fS_continent_bmtt,
                    Sstd2_lognorm_fS_cbmtt1)
                    # Sstd2_lognorm_fS_clim_cbmtt0,
                    # Sstd2_lognorm_fS_clim_cbmtt1) 

# save(Sstd2_lognorm_fragSize,
#      Sstd2_lognorm_fS_clim_cont,
#      Sstd2_lognorm_fS_clim_cont_bmtt,
#      Sstd2_lognorm_fS_cbmtt0,
#      Sstd2_lognorm_fS_continent0,
#      Sstd2_lognorm_fS_continent1,
#      Sstd2_lognorm_fS_continent_bmtt,
#      Sstd2_lognorm_fS_cbmtt1,
#      Sstd2_lognorm_fS_clim_cbmtt0,
#      Sstd2_lognorm_fS_clim_cbmtt1,
#      file = '~/Dropbox/1current/fragmentation_synthesis/results/Sstd2_model_explore.Rdata')
Sstd2_lognorm_fS_cont2 <- update(Sstd2_lognorm_fragSize,
                                  formula. = ~ 1 + (c.lfs*climate | cbmtt / dataset_label), 
                                  newdata = frag,
                                  # prior = rp,
                                  family = 'lognormal',
                                  cores = 4, chains = 4)

Sn_lognorm_fragSize <- brm(S_n ~ c.lfs + (c.lfs | dataset_label), 
                           data = frag,
                           prior = rp,
                           family = 'lognormal',
                           cores = 4, chains = 4)

Sn_lognorm_fragSize2 <- update(Sn_lognorm_fragSize,
                                  formula. = ~ c.lfs*cbmtt + (c.lfs | dataset_label),
                                  newdata = frag,
                                  # prior = rp,
                                  family = 'lognormal',
                                  cores = 4, chains = 4)

Sn_lognorm_fragSize4 <- update(Sn_lognorm_fragSize,
                               formula. = ~ c.lfs + (c.lfs | cbmtt/dataset_label), 
                               newdata = frag,
                               # prior = rp,
                               family = 'lognormal',
                               cores = 4, chains = 4)



Scov_lognorm_fragSize <- brm(S_cov ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag,
                             prior = rp,
                             family = hurdle_lognormal(),
                             cores = 4, chains = 4)

Scov_lognorm_fragSize2 <- update(Scov_lognorm_fragSize,
                                 formula. = ~ c.lfs*cbmtt + (c.lfs | dataset_label),
                                 newdata = frag,
                                 # prior = rp,
                                 family = hurdle_lognormal(),
                                 cores = 4, chains = 4)

Scov_lognorm_fragSize4 <- update(Scov_lognorm_fragSize,
                               formula. = ~ c.lfs + (c.lfs | cbmtt/dataset_label), 
                               newdata = frag,
                               # prior = rp,
                               family = hurdle_lognormal(),
                               cores = 4, chains = 4)

S_PIE_lognorm_fragSize <- brm(S_PIE ~ c.lfs + (c.lfs | dataset_label), 
                              data = frag,
                              prior = rp,
                              family = 'lognormal',
                              cores = 4, chains = 4)

# S_PIE_lognorm_fragSize3 <- update(S_PIE_lognorm_fragSize,
#                                  formula. = ~ c.lfs + (c.lfs | bmtt/dataset_label), 
#                                  newdata = frag,
#                                  # prior = rp,
#                                  family = 'lognormal',
#                                  cores = 4, chains = 4)

S_PIE_lognorm_fragSize4 <- update(S_PIE_lognorm_fragSize,
                                  formula. = ~ c.lfs + (c.lfs | cbmtt/dataset_label), 
                                  newdata = frag,
                                  # prior = rp,
                                  family = 'lognormal',
                                  cores = 4, chains = 4)



S_chao_lognorm_fragSize <- brm(S_chao ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag,
                               prior = rp,
                               family = 'lognormal',
                               cores = 4, chains = 4)

# S_chao_lognorm_fragSize3 <- update(S_chao_lognorm_fragSize,
#                                    formula. = ~ c.lfs + (c.lfs | bmtt/dataset_label), 
#                                    newdata = frag,
#                                    # prior = rp,
#                                    family = 'lognormal',
#                                    cores = 4, chains = 4)

S_chao_lognorm_fragSize4 <- update(S_chao_lognorm_fragSize,
                                  formula. = ~ c.lfs + (c.lfs | cbmtt/dataset_label), 
                                  newdata = frag,
                                  # prior = rp,
                                  family = 'lognormal',
                                  cores = 4, chains = 4)



Nstd_lognorm_fragSize <- brm(N_std ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag,
                             prior = rp,
                             family = 'lognormal',
                             cores = 4, chains = 4)

# Nstd_lognorm_fragSize3 <- update(Nstd_lognorm_fragSize,
#                                  formula. = ~ c.lfs + (c.lfs | bmtt/dataset_label), 
#                                  newdata = frag,
#                                  # prior = rp,
#                                  family = 'lognormal',
#                                  cores = 4, chains = 4)

Nstd_lognorm_fragSize4 <- update(Nstd_lognorm_fragSize,
                                   formula. = ~ c.lfs + (c.lfs | cbmtt/dataset_label), 
                                   newdata = frag,
                                   # prior = rp,
                                   family = 'lognormal',
                                   cores = 4, chains = 4)



save(Sstd1_lognorm_fragSize, 
     Sstd2_lognorm_fragSize,
     Sstd2_lognorm_fragSize4,
     Sn_lognorm_fragSize,
     Sn_lognorm_fragSize4,
     S_PIE_lognorm_fragSize,
     S_PIE_lognorm_fragSize4,
     Scov_lognorm_fragSize,
     Scov_lognorm_fragSize4,
     S_chao_lognorm_fragSize,
     S_chao_lognorm_fragSize4,
     Nstd_lognorm_fragSize,
     Nstd_lognorm_fragSize4,
     file = '~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms.Rdata')

waic(Sstd2_lognorm_fragSize,
     Sstd2_lognorm_fS_cont1,
     Sstd2_lognorm_fragSize4) # slight improvement (3.7) for cbmtt as grouping variable

pp_check(Sstd2_lognorm_fragSize4) +
  scale_x_continuous(trans = 'log', 
                     breaks = c(1,8,16,32,64,128))

waic(Sn_lognorm_fragSize,
     Sn_lognorm_fragSize3,
     Sn_lognorm_fragSize4) # including group slightly worse 2.5
# 
pp_check(Sn_lognorm_fragSize4) +
  scale_x_continuous(trans = 'log',
                     breaks = c(1,8,16,32,64,128))

waic(S_chao_lognorm_fragSize,
     S_chao_lognorm_fragSize4) # slightly worse 3.3
# 
pp_check(S_chao_lognorm_fragSize4) +
  scale_x_continuous(trans = 'log',
                     breaks = c(1,8,16,32,64,128))

waic(Scov_lognorm_fragSize,
     Scov_lognorm_fragSize4) # slightly worse 2.9
# 
pp_check(Scov_lognorm_fragSize) +
  scale_x_continuous(trans = 'log',
                     breaks = c(1,8,16,32,64,128))

waic(S_PIE_lognorm_fragSize,
     S_PIE_lognorm_fragSize4) # no change 0.01
# 
pp_check(S_PIE_lognorm_fragSize) +
  scale_x_continuous(trans = 'log',
                     breaks = c(1,8,16,32,64,128))

waic(Nstd_lognorm_fragSize,
     Nstd_lognorm_fragSize4) # no change 0.3

pp_check(Nstd_lognorm_fragSize4) +
  scale_x_continuous(trans = 'log') # 