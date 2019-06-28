rm(list=ls())
# code to fit refit models for fragmentation synthesis to examine sensitivity to decisions made
# in the calculation (standardisation) of the metrics

library(tidyverse)
library(brms)

# list of files
files = list.files(path = '~/Dropbox/Frag Database (new)/files_datapaper/Analysis/Sensitivity_analysis/priority/',
                   pattern = 'biodiv_frag_fcont')

for(i in 1:length(files)){
  # get the data
  file_2_get = paste0('~/Dropbox/Frag Database (new)/files_datapaper/Analysis/Sensitivity_analysis/priority/',
                     files[i])
  frag <- read_csv(file_2_get)  
  
  # add mean-centred log(fragment.size)
  frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))
  
  # set some weakly regularising priors (as used for reference analysis)
  rp <- c(prior(normal(2,2), class = Intercept),
          prior(normal(0,1), class = b),
          prior(exponential(1), class = sd))  
  
  Sstd2_lognorm_fragSize <- brm(S_std_2 ~ c.lfs + (c.lfs | dataset_label), 
                                # fit to data with variation in frag_size_num
                                data = frag,
                                prior = rp,
                                family = 'lognormal', # our standardised richness are not integer values
                                cores = 4, chains = 4)

  Sn_lognorm_fragSize <- brm(S_n ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag,
                             prior = rp,
                             family = 'lognormal',
                             cores = 4, chains = 4)
  
  Scov_lognorm_fragSize <- brm(S_cov ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag,
                               prior = rp,
                               family = hurdle_lognormal(),
                               cores = 4, chains = 4)

  S_PIE_lognorm_fragSize <- brm(S_PIE ~ c.lfs + (c.lfs | dataset_label), 
                                data = frag,
                                prior = rp,
                                family = 'lognormal',
                                cores = 4, chains = 4)
  S_chao_lognorm_fragSize <- brm(S_chao ~ c.lfs + (c.lfs | dataset_label), 
                                 data = frag,
                                 prior = rp,
                                 family = 'lognormal',
                                 cores = 4, chains = 4)
  
  Nstd_lognorm_fragSize <- brm(N_std ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag,
                               prior = rp,
                               family = 'lognormal',
                               cores = 4, chains = 4)
  name_2_save = paste0('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_', 
                       strsplit(files[i], split = '.csv')[[1]])
  
  save(#Sstd1_lognorm_fragSize, 
    Sstd2_lognorm_fragSize,
    Sn_lognorm_fragSize,
    S_PIE_lognorm_fragSize,
    Scov_lognorm_fragSize,
    S_chao_lognorm_fragSize,
    Nstd_lognorm_fragSize,
    file = name_2_save)
}
