# set paths and load libraries using 0_init_dirs_load_packages.R

# code to fit refit models for fragmentation synthesis to examine sensitivity to decisions made
# in the imputation of habitat area for "continuous habitat" and in the treatment of non-integer abundance values

# list of files
files = list.files(path = paste0(path2wd, 'intermediate_results/'),
                   pattern = 'biodiv_frag_fcont')

# reference case are models already fit to 2_biodiv_frag_fcont_10_mabund_as_is.csv
files = files[-which(files=='2_biodiv_frag_fcont_10_mabund_as_is.csv')]

for(i in 2:length(files)){
  # get the data
  file_2_get = paste0(path2wd, 'intermediate_results/',
                     files[i])
  frag <- read_csv(file_2_get)  
  
  # add mean-centred log(fragment.size)
  frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))
  

  Sstd_lognorm_fragSize <- brm(S_std_mean ~ c.lfs + (c.lfs | dataset_label), 
                                # fit to data with variation in frag_size_num
                                data = frag %>% filter(S_std_mean > 0),
                                family = 'lognormal', # our standardised richness are not integer values
                                cores = 4, chains = 4)

  S_PIE_lognorm_fragSize <- brm(S_PIE_mean ~ c.lfs + (c.lfs | dataset_label), 
                                data = frag %>% filter(S_PIE_mean > 0),
                                family = 'lognormal',
                                cores = 4, chains = 4)

  Nstd_lognorm_fragSize <- brm(N_std ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag,
                               family = 'lognormal',
                               cores = 4, chains = 4)
  name_2_save = paste0(path2wd, 'extended_data_figs_tabs/', 
                       strsplit(files[i], split = '.csv')[[1]], '_modelFit.Rdata')
  
  save(Sstd_lognorm_fragSize,
    S_PIE_lognorm_fragSize,
    Nstd_lognorm_fragSize,
    file = name_2_save)
}
