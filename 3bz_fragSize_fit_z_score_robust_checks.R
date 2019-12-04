# set paths and load libraries using 0_init_dirs_load_packages.R

# code to fit refit models for fragmentation synthesis to examine sensitivity to decisions made
# in the calculation (standardisation) of the metrics

# list of files
files = list.files(path = path2data,
                   pattern = 'biodiv_frag_fcont')
# I have already fit modesl to 2_biodiv_frag_fcont_10_mabund_as_is.csv
files = files[-which(files=='2_biodiv_frag_fcont_10_mabund_as_is.csv')]

for(i in 2:length(files)){
  # get the data
  file_2_get = paste0(path2data,
                      files[i])
  frag <- read_csv(file_2_get)  
  
  # add mean-centred log(fragment.size)
  frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))
  
  #reference analysis used default priors
  # rp <- c(prior(normal(2,2), class = Intercept),
  #         prior(normal(0,1), class = b),
  #         prior(exponential(1), class = sd))  
  
  z_Sstd_studT_fragSize <- brm(z_S_std ~ c.lfs + (c.lfs | dataset_label), 
                              # some z-scores are infinite due to sd(expected) = 0
                              data = frag %>% filter(S_std_mean>0 & 
                                                       !is.na(z_S_std) & 
                                                       !is.infinite(z_S_std)),
                              # prior = rp,
                              # spike in values near zero is not well described with gaussian error
                              # alternate use student T
                              family = student(),
                              cores = 4, 
                              chains = 4
  )
  
  z_Sn_studT_fragSize <- brm(z_S_n ~ c.lfs + (c.lfs | dataset_label), 
                             # some z-scores are infinite due to sd(expected) = 0
                             data = frag %>% filter(S_n_mean>0 & 
                                                      !is.na(z_S_n) & 
                                                      !is.infinite(z_S_n)),
                             # prior = rp,
                             # spike in values near zero is not well described with gaussian error
                             # alternate use student T
                             family = student(),
                             cores = 4, 
                             chains = 4
  )
  
  
  z_Scov_studT_fragSize <- brm(z_S_cov ~ c.lfs + (c.lfs | dataset_label), 
                               # some z-scores are infinite due to sd(expected) = 0
                               data = frag %>% filter(S_cov_mean>0 & 
                                                        !is.na(z_S_cov) & 
                                                        !is.infinite(z_S_cov)),
                               # prior = rp,
                               # spike in values near zero is not well described with gaussian error
                               # alternate use student T
                               family = student(),
                               cores = 4, 
                               chains = 4
  )
  
  z_S_PIE_studT_fragSize <- brm(z_S_PIE ~ c.lfs + (c.lfs | dataset_label), 
                                # some z-scores are infinite due to sd(expected) = 0
                                data = frag %>% filter(S_PIE_mean>0 & 
                                                         !is.na(z_S_PIE) & 
                                                         !is.infinite(z_S_PIE)),
                                # prior = rp,
                                # spike in values near zero is not well described with gaussian error
                                # alternate use student T
                                family = student(),
                                cores = 4, 
                                chains = 4
  )
  
  z_S_chao_studT_fragSize <- brm(z_S_chao ~ c.lfs + (c.lfs | dataset_label), 
                                 # some z-scores are infinite due to sd(expected) = 0
                                 data = frag %>% filter(S_chao_mean>0 & 
                                                          !is.na(z_S_chao) & 
                                                          !is.infinite(z_S_chao)),
                                 # prior = rp,
                                 # what about Student's distribution?
                                 family = student(),
                                 cores = 4, 
                                 chains = 4
  )
  name_2_save = paste0('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_', 
                       strsplit(files[i], split = '.csv')[[1]], '.Rdata')
  
  save(#Sstd1_lognorm_fragSize, 
    z_Sstd_studT_fragSize,
    z_Sn_studT_fragSize,
    z_Scov_studT_fragSize,
    z_S_PIE_studT_fragSize,
    z_S_chao_studT_fragSize,
    file = name_2_save)
}
