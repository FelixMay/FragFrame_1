# code to fit models for fragmentation synthesis
# so far: bayesian framework for approximately ML-like results (i.e.,
# with non-informative priors - against the science of Gelman)
library(tidyverse)
library(brms)

# load the data
frag <- read_csv('~/Dropbox/Frag Database (new)/files_datapaper/Analysis/2_biodiv_frag_fcont_10_mabund_as_is.csv')

# visual inspection 1
frag %>% 
  ggplot() +
  geom_histogram(aes(x = frag_size_num)) +
  scale_x_continuous(trans = 'log')

# visual inspection 2: missing some S_std_1 and S_std_2 values
frag %>% 
  ggplot() +
  geom_point(aes(x = frag_size_num, y = S_std_1)) +
  stat_smooth(method = 'lm', se = F,
              aes(x = frag_size_num, y = S_std_1),
              col = 'black') +
  geom_point(aes(x = frag_size_num, y = S_std_2),
             col = 'red') +
  stat_smooth(method = 'lm', se = F,
              aes(x = frag_size_num, y = S_std_2),
              col = 'red') +
  scale_x_continuous(trans = 'log') +
  scale_y_continuous(trans = 'log') +
  theme_bw() +
  theme(legend.position = 'none')


##--create some covariates for easier workflow--
# mean-centred log(fragment.size)
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))
# add log-transformed (to fit multivariate normal?)
frag$lSstd1 <- log(frag$S_std_1)
frag$lSstd2 <- log(frag$S_std_2)
frag$lSn <- log(frag$S_n)
frag$lScov <- log(frag$S_cov)
frag$lS_PIE <- log(frag$S_PIE)
frag$lS_chao <- log(frag$S_chao)
frag$lNstd <- log(frag$N_std)

# how steep is the slope? Think about regularising prior (plot for completeness)
with(frag, summary(lm(lSstd1~c.lfs)))

#----- simplest model: diversity as a function of fragment size; allow fragment size to vary by study----
get_prior(S_std_2 ~ c.lfs + (c.lfs | dataset_label),
          family = lognormal(),
          data = frag)
# set some weakly regularising priors
rp <- c(prior(normal(0,2), class = Intercept),
        prior(normal(0,1), class = b),
        prior(exponential(1), class = sd))

# not sure which of the standardised richness values to model, do both
Sstd1_lognorm_fragSize <- brm(S_std_1 ~ c.lfs + (c.lfs | dataset_label), 
                              # fit to data with variation in frag_size_num
                              data = frag,
                              prior = rp,
                              family = 'lognormal',
                              cores = 4, chains = 4)

pp_check(Sstd1_lognorm_fragSize) +
  scale_x_continuous(trans = 'log', 
                     breaks = c(1,8,16,32,64,128))

Sstd2_lognorm_fragSize <- brm(S_std_2 ~ c.lfs + (c.lfs | dataset_label), 
                              # fit to data with variation in frag_size_num
                              data = frag,
                              prior = rp,
                              family = 'lognormal',
                              cores = 4, chains = 4)

pp_check(Sstd2_lognorm_fragSize) +
  scale_x_continuous(trans = 'log', 
                     breaks = c(1,8,16,32,64,128))

Sn_lognorm_fragSize <- brm(S_n ~ c.lfs + (c.lfs | dataset_label), 
                           data = frag,
                           prior = rp,
                           family = 'lognormal',
                           cores = 4, chains = 4)

pp_check(Sn_lognorm_fragSize) +
  scale_x_continuous(trans = 'log', 
                     breaks = c(1,8,16,32,64,128))

Scov_lognorm_fragSize <- brm(S_cov ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag,
                             prior = rp,
                             family = hurdle_lognormal(),
                             cores = 4, chains = 4)

pp_check(Scov_lognorm_fragSize) +
  scale_x_continuous(trans = 'log', 
                     breaks = c(1,8,16,32,64,128))

S_PIE_lognorm_fragSize <- brm(S_PIE ~ c.lfs + (c.lfs | dataset_label), 
                              data = frag,
                              prior = rp,
                              family = 'lognormal',
                              cores = 4, chains = 4)

pp_check(S_PIE_lognorm_fragSize) +
  scale_x_continuous(trans = 'log', 
                     breaks = c(1,8,16,32,64,128))

S_chao_lognorm_fragSize <- brm(S_chao ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag,
                               prior = rp,
                               family = 'lognormal',
                               cores = 4, chains = 4)

pp_check(S_chao_lognorm_fragSize) +
  scale_x_continuous(trans = 'log', 
                     breaks = c(1,8,16,32,64,128))

# note sure which abundance to model: do both?
# N_lognorm_fragSize <- brm(N ~ c.lfs + (c.lfs | dataset_label), 
#                                data = frag,
#                                prior = rp,
#                                family = 'lognormal',
#                                cores = 4, chains = 4)

Nstd_lognorm_fragSize <- brm(N_std ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag,
                             prior = rp,
                             family = 'lognormal',
                             cores = 4, chains = 4)

pp_check(Nstd_lognorm_fragSize) +
  scale_x_continuous(trans = 'log') # 

save(Sstd1_lognorm_fragSize, 
     Sstd2_lognorm_fragSize,
     Sn_lognorm_fragSize,
     S_PIE_lognorm_fragSize,
     Scov_lognorm_fragSize,
     S_chao_lognorm_fragSize,
     Nstd_lognorm_fragSize,
     file = '~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms.Rdata')




