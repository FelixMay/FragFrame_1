# set paths and load libraries using 0_init_dirs_load_packages.R

# code to fit refit models for fragmentation synthesis to examine sensitivity to decisions made
# in the calculation (standardisation) of the metrics

# list of files
files = list.files(path = paste0(path2wd, 'intermediate_results/'),
                   pattern = 'biodiv_frag_fcont')

# reference case are models already fit to 2_biodiv_frag_fcont_10_mabund_as_is.csv
files = files[-which(files=='2_biodiv_frag_fcont_10_mabund_as_is.csv')]

for(i in 1:length(files)){
  # get the data
  file_2_get = paste0(path2wd, 'intermediate_results/',
                     files[i])
  frag <- read_csv(file_2_get)  
  
  # add mean-centred log(fragment.size)
  frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))
  
  #reference analysis used default priors
  # rp <- c(prior(normal(2,2), class = Intercept),
  #         prior(normal(0,1), class = b),
  #         prior(exponential(1), class = sd))  
  
  Sstd2_lognorm_fragSize <- brm(S_std2_mean ~ c.lfs + (c.lfs | dataset_label), 
                                # fit to data with variation in frag_size_num
                                data = frag %>% filter(S_std1_mean > 0),
                                # prior = rp,
                                family = 'lognormal', # our standardised richness are not integer values
                                cores = 4, chains = 4)

  Sn_lognorm_fragSize <- brm(S_n_mean ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag %>% filter(S_n_mean > 0),
                             # prior = rp,
                             family = 'lognormal',
                             cores = 4, chains = 4)
  
  Scov_lognorm_fragSize <- brm(S_cov_mean ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag %>% filter(S_cov_mean > 0),
                               # prior = rp,
                               family = hurdle_lognormal(),
                               cores = 4, chains = 4)

  S_PIE_lognorm_fragSize <- brm(S_PIE_mean ~ c.lfs + (c.lfs | dataset_label), 
                                data = frag %>% filter(S_PIE_mean > 0),
                                # prior = rp,
                                family = 'lognormal',
                                cores = 4, chains = 4)
  S_chao_lognorm_fragSize <- brm(S_chao_mean ~ c.lfs + (c.lfs | dataset_label), 
                                 data = frag %>% filter(S_chao_mean > 0),
                                 # prior = rp,
                                 family = 'lognormal',
                                 cores = 4, chains = 4)
  
  Nstd_lognorm_fragSize <- brm(N_std ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag,
                               # prior = rp,
                               family = 'lognormal',
                               cores = 4, chains = 4)
  name_2_save = paste0(path2wd, 'intermediate_results/', 
                       strsplit(files[i], split = '.csv')[[1]], '_modelFits.Rdata')
  
  save(#Sstd1_lognorm_fragSize, 
    Sstd2_lognorm_fragSize,
    Sn_lognorm_fragSize,
    S_PIE_lognorm_fragSize,
    Scov_lognorm_fragSize,
    S_chao_lognorm_fragSize,
    Nstd_lognorm_fragSize,
    file = name_2_save)
}
