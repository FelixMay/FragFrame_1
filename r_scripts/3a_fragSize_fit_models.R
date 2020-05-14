# need to execute 0_init_dirs_load_packages.R first

# code to fit models for fragmentation synthesis
# for this analysis I am sticking with the default priors for simplicity 
# they are weakly regularising

# load the data
frag <- read_csv(paste0(path2wd, 'intermediate_results/2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))


#----- simplest model: diversity as a function of fragment size; 
# allow fragment size (slope) to vary by study (varying intercept)----

# results presented in ms are fit with default, weakly regularising priors
# that are informed by the data (see brms documentation)

# sample effort standardised species richness 
Sstd_lognorm_fragSize <- brm(S_std_mean ~ c.lfs + (c.lfs | dataset_label), 
                              # two observations have zero for the response, remove for lognormal distribution
                              data = frag %>% filter(S_std_mean>0),
                              family = lognormal(), # our standardised richness are not integer values
                              cores = 4, chains = 4,
                             iter = 4000, thin = 2)

# refit without the studies with pooled sampling designs
Sstd_lognorm_fragSize_pool <- brm(S_std_mean ~ c.lfs + (c.lfs | dataset_label), 
                              # remove pooled designs
                              data = frag %>% filter(sample_design!='pooled'),
                              family = lognormal(), # our standardised richness are not integer values
                              cores = 4, chains = 4,
                              iter = 4000, thin = 2,
                              control = list(adapt_delta = 0.95))

# richness standardised by number of individuals
Sn_lognorm_fragSize <- brm(S_n_mean ~ c.lfs + (c.lfs | dataset_label), 
                           data = frag %>% filter(S_n_mean > 0), # for lognormal distribution
                           family = lognormal(),
                           cores = 4, chains = 4,
                           iter = 4000, thin = 2)

# refit without the studies with pooled sampling designs
Sn_lognorm_fragSize_pool <- brm(S_n_mean ~ c.lfs + (c.lfs | dataset_label), 
                           data = frag %>% filter(S_n_mean > 0 & sample_design!='pooled'),
                           family = lognormal(),
                           cores = 4, chains = 4,
                           iter = 4000, thin = 2)

# coverage standardised richness
Scov_lognorm_fragSize <- brm(S_cov_mean ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag %>% filter(S_cov_mean > 0),
                             family = lognormal(),
                             cores = 4, chains = 4,
                             iter = 4000, thin = 2)

# refit without studies with pooled sampling designs
Scov_lognorm_fragSize_pool <- brm(S_cov_mean ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag %>% filter(S_cov_mean > 0 & sample_design!='pooled'),
                             family = lognormal(),
                             cores = 4, chains = 4,
                             iter = 4000, thin = 2)

# effective number of species conversion of the probability of interspecific encounter
S_PIE_lognorm_fragSize <- brm(S_PIE_mean ~ c.lfs + (c.lfs | dataset_label), 
                              data = frag %>% filter(!is.na(S_PIE_mean)),
                              family = lognormal(),
                              cores = 4, chains = 4,
                              iter = 4000, thin = 2)

# refit without studies with pooled sampling designs
S_PIE_lognorm_fragSize_pool <- brm(S_PIE_mean ~ c.lfs + (c.lfs | dataset_label), 
                              data = frag %>% filter(!is.na(S_PIE_mean) & sample_design!='pooled'),
                              family = lognormal(),
                              cores = 4, chains = 4,
                              iter = 4000, thin = 2)

S_chao_lognorm_fragSize <- brm(S_chao_mean ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag %>% filter(S_chao_mean > 0),
                               family = lognormal(),
                               cores = 4, chains = 4,
                               iter = 4000, thin = 2)

# refit without the studies with pooled sampling designs
S_chao_lognorm_fragSize_pool <- brm(S_chao_mean ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag %>% filter(S_chao_mean > 0 & sample_design!='pooled'),
                               family = lognormal(),
                               cores = 4, chains = 4,
                               iter = 4000, thin = 2)


Nstd_lognorm_fragSize <- brm(N_std ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag,
                             family = lognormal(),
                             cores = 4, chains = 4,
                             iter = 4000, thin = 2)


Nstd_lognorm_fragSize_pool <- brm(N_std ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag %>% filter(sample_design!='pooled'),
                             family = lognormal(),
                             cores = 4, chains = 4,
                             iter = 4000, thin = 2)


save(Sstd_lognorm_fragSize,
     Sn_lognorm_fragSize,
     S_PIE_lognorm_fragSize,
     Scov_lognorm_fragSize,
     S_chao_lognorm_fragSize,
     Nstd_lognorm_fragSize,
     file = paste0(path2wd, 'main_results/fragSize_ref.Rdata'))

save(Sstd_lognorm_fragSize_pool,
     Sn_lognorm_fragSize_pool,
     S_PIE_lognorm_fragSize_pool,
     Scov_lognorm_fragSize_pool,
     S_chao_lognorm_fragSize_pool,
     Nstd_lognorm_fragSize_pool,
     file = paste0(path2wd, 'main_results/fragSize_ref_pool.Rdata'))