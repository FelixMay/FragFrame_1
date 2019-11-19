# need to execute 0_init_dirs_load_packages.R first

# code to fit models for fragmentation synthesis
# for this analysis I am sticking with the default priors for simplicity 
# they are weakly regularising

# load the data
frag <- read_csv(paste0(path2data, '2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))


#----- simplest model: diversity as a function of fragment size; allow fragment size to vary by study----

# results presented in ms are fit with default priors

# alternate to default priors: weakly regularising priors
# rp <- c(prior(normal(0, 1), class = Intercept),
#         prior(normal(0, 0.5), class = b),
#         prior(exponential(1), class = sd))
 

# sample effort standardised species richness 
Sstd_lognorm_fragSize <- brm(S_std_mean ~ c.lfs + (c.lfs | dataset_label), 
                              # fit to data with variation in frag_size_num
                              data = frag %>% filter(S_std_mean>0),
                              #prior = rp,
                              family = 'lognormal', # our standardised richness are not integer values
                              cores = 4, chains = 4)

# refit without the studies with pooled sampling designs
Sstd_lognorm_fragSize_pool <- brm(S_std_mean ~ c.lfs + (c.lfs | dataset_label), 
                              # fit to data with variation in frag_size_num
                              data = frag %>% filter(sample_design!='pooled'),
                              #prior = rp,
                              family = 'lognormal', # our standardised richness are not integer values
                              cores = 4, chains = 4,
                              control = list(adapt_delta = 0.95))

# richness standardised by number of individuals
Sn_lognorm_fragSize <- brm(S_n_mean ~ c.lfs + (c.lfs | dataset_label), 
                           data = frag %>% filter(S_n_mean > 0),
                           #prior = rp,
                           family = 'lognormal',
                           cores = 4, chains = 4)

# refit without the studies with pooled sampling designs
Sn_lognorm_fragSize_pool <- brm(S_n_mean ~ c.lfs + (c.lfs | dataset_label), 
                           data = frag %>% filter(S_n_mean > 0 & sample_design!='pooled'),
                           #prior = rp,
                           family = 'lognormal',
                           cores = 4, chains = 4)

# coverage standardised richness
Scov_lognorm_fragSize <- brm(S_cov_mean ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag %>% filter(S_cov_mean > 0),
                             #prior = rp,
                             family = lognormal(),
                             cores = 4, chains = 4)

# refit without studies with pooled sampling designs
Scov_lognorm_fragSize_pool <- brm(S_cov_mean ~ c.lfs + (c.lfs | dataset_label), 
                             data = frag %>% filter(S_cov_mean > 0 & sample_design!='pooled'),
                             #prior = rp,
                             family = lognormal(),
                             cores = 4, chains = 4)

# effective number of species conversion of the probability of interspecific encounter
S_PIE_lognorm_fragSize <- brm(S_PIE_mean ~ c.lfs + (c.lfs | dataset_label), 
                              data = frag %>% filter(!is.na(S_PIE_mean)),
                              #prior = rp,
                              family = 'lognormal',
                              cores = 4, chains = 4)

# refit without studies with pooled sampling designs
S_PIE_lognorm_fragSize_pool <- brm(S_PIE_mean ~ c.lfs + (c.lfs | dataset_label), 
                              data = frag %>% filter(!is.na(S_PIE_mean) & sample_design!='pooled'),
                              #prior = rp,
                              family = 'lognormal',
                              cores = 4, chains = 4)

S_chao_lognorm_fragSize <- brm(S_chao_mean ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag %>% filter(S_chao_mean > 0),
                               #prior = rp,
                               family = 'lognormal',
                               cores = 4, chains = 4)

# refit without the studies with pooled sampling designs
S_chao_lognorm_fragSize_pool <- brm(S_chao_mean ~ c.lfs + (c.lfs | dataset_label), 
                               data = frag %>% filter(S_chao_mean > 0 & sample_design!='pooled'),
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


save(Sstd_lognorm_fragSize,
     Sstd_lognorm_fragSize_pool,
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
     file = '~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms_ref_revision.Rdata')
