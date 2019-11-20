# need to execute 0_init_dirs_load_packages.R first

# code to fit models with 2-way interaction for fragmentation synthesis
# so far: bayesian framework for approximately ML-like results (i.e.,
# with non-informative priors - against the science of Gelman)

frag <- read_csv(paste0(path2data, '2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

meta <- read.csv(paste0(path2meta, 'new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')

##--create some covariates for easier workflow--

# # set the reference levels for the categorical covariates of interest
frag %>% distinct(Matrix.category)
frag$Matrix.category <- factor(frag$Matrix.category, levels = c('light filter', 'intermediate', 'harsh filter'))
frag %>% distinct(time.since.fragmentation)
frag$time.since.fragmentation <- factor(frag$time.since.fragmentation,
                                         levels = c("Recent (less than 20 years)", "Intermediate (20-100 years)", "long (100+ years)"))

frag %>% distinct(continent8)

# load fragSize only models
load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms_ref_revision.Rdata')
# two-way interactions: matrix permeability first
Sstd_ln_fS_matrix <- update(Sstd_lognorm_fragSize,
                            formula. = ~ c.lfs * Matrix.category + (c.lfs | dataset_label),
                            newdata = frag %>% filter(S_std_mean > 0), 
                            cores = 4)

Sn_ln_fS_matrix <- update(Sn_lognorm_fragSize, 
                          formula. = ~ c.lfs * Matrix.category + (c.lfs | dataset_label),
                          newdata = frag %>% filter(S_n_mean > 0),
                          cores = 4)

Scov_ln_fS_matrix <- update(Scov_lognorm_fragSize, 
                          formula. = ~ c.lfs * Matrix.category + (c.lfs | dataset_label),
                          newdata = frag %>% filter(S_cov_mean > 0), 
                          cores = 4)

Schao_ln_fS_matrix <- update(S_chao_lognorm_fragSize, 
                             formula. = ~ c.lfs * Matrix.category + (c.lfs | dataset_label),
                             newdata = frag %>% filter(S_chao_mean > 0), 
                             cores = 4)

S_PIE_ln_fS_matrix <- update(S_PIE_lognorm_fragSize, 
                             formula. = ~ c.lfs * Matrix.category + (c.lfs | dataset_label),
                             newdata = frag %>% filter(S_PIE_mean > 0),
                             cores = 4)

N_std_ln_fS_matrix <- update(Nstd_lognorm_fragSize, 
                             formula. = ~ c.lfs * Matrix.category + (c.lfs | dataset_label),
                             newdata = frag, 
                             cores = 4)  

# repeat for taxa
Sstd_ln_fS_taxa <- update(Sstd_lognorm_fragSize, 
                          formula. = ~ c.lfs * taxa + (c.lfs | dataset_label),
                          newdata = frag %>% filter(S_std_mean > 0), 
                          cores = 4)

Sn_ln_fS_taxa <- update(Sn_lognorm_fragSize, 
                        formula. = ~ c.lfs * taxa + (c.lfs | dataset_label),
                        newdata = frag %>% filter(S_n_mean > 0), 
                        cores = 4)

Scov_ln_fS_taxa <- update(Scov_lognorm_fragSize, 
                          formula. = ~ c.lfs * taxa + (c.lfs | dataset_label),
                          newdata = frag %>% filter(S_cov_mean > 0), 
                          cores = 4)

Schao_ln_fS_taxa <- update(S_chao_lognorm_fragSize, 
                           formula. = ~ c.lfs * taxa + (c.lfs | dataset_label),
                           newdata = frag %>% filter(S_chao_mean > 0), 
                           cores = 4)

S_PIE_ln_fS_taxa <- update(S_PIE_lognorm_fragSize, 
                           formula. = ~ c.lfs * taxa + (c.lfs | dataset_label),
                           newdata = frag %>% filter(S_PIE_mean > 0), 
                           cores = 4)

N_std_ln_fS_taxa <- update(Nstd_lognorm_fragSize, 
                           formula. = ~ c.lfs * taxa + (c.lfs | dataset_label),
                           newdata = frag, 
                           cores = 4)  

# repeat for time since fragmentation
Sstd_ln_fS_tsf <- update(Sstd_lognorm_fragSize, 
                         formula. = ~ c.lfs * time.since.fragmentation + (c.lfs | dataset_label),
                         newdata = frag %>% filter(S_std_mean > 0), 
                         cores = 4)

Sn_ln_fS_tsf <- update(Sn_lognorm_fragSize, 
                       formula. = ~ c.lfs * time.since.fragmentation + (c.lfs | dataset_label),
                       newdata = frag %>% filter(S_n_mean > 0), 
                       cores = 4)

Scov_ln_fS_tsf <- update(Scov_lognorm_fragSize, 
                         formula. = ~ c.lfs * time.since.fragmentation + (c.lfs | dataset_label),
                         newdata = frag %>% filter(S_cov_mean > 0),
                         cores = 4)

Schao_ln_fS_tsf <- update(S_chao_lognorm_fragSize, 
                          formula. = ~ c.lfs * time.since.fragmentation + (c.lfs | dataset_label),
                          newdata = frag %>% filter(S_chao_mean > 0), 
                          cores = 4)

S_PIE_ln_fS_tsf <- update(S_PIE_lognorm_fragSize, 
                          formula. = ~ c.lfs * time.since.fragmentation + (c.lfs | dataset_label),
                          newdata = frag %>% filter(S_PIE_mean > 0), 
                          cores = 4)

N_std_ln_fS_tsf <- update(Nstd_lognorm_fragSize, 
                          formula. = ~ c.lfs * time.since.fragmentation + (c.lfs | dataset_label),
                          newdata = frag, 
                          cores = 4)  

# repeat for biome
Sstd_ln_fS_biome <- update(Sstd_lognorm_fragSize, 
                          formula. = ~ c.lfs * biome + (c.lfs | dataset_label),
                          newdata = frag %>% filter(S_std_mean > 0),
                          cores = 4)

Sn_ln_fS_biome <- update(Sn_lognorm_fragSize, 
                       formula. = ~ c.lfs * biome + (c.lfs | dataset_label),
                       newdata = frag %>% filter(S_n_mean > 0),
                       cores = 4)

Scov_ln_fS_biome <- update(Scov_lognorm_fragSize, 
                         formula. = ~ c.lfs * biome + (c.lfs | dataset_label),
                         newdata = frag %>% filter(S_cov_mean > 0),
                         cores = 4)

Schao_ln_fS_biome <- update(S_chao_lognorm_fragSize, 
                          formula. = ~ c.lfs * biome + (c.lfs | dataset_label),
                          newdata = frag %>% filter(S_chao_mean > 0),
                          cores = 4)

S_PIE_ln_fS_biome <- update(S_PIE_lognorm_fragSize, 
                          formula. = ~ c.lfs * biome + (c.lfs | dataset_label),
                          newdata = frag %>% filter(S_PIE_mean > 0),
                          cores = 4)

N_std_ln_fS_biome <- update(Nstd_lognorm_fragSize, 
                            formula. = ~ c.lfs * biome + (c.lfs | dataset_label),
                            newdata = frag,
                            cores = 4)
# region
Sstd_ln_fS_region <- update(Sstd_lognorm_fragSize,
                            formula. = ~ c.lfs * continent8 + (c.lfs | dataset_label),
                            newdata = frag %>% filter(S_std_mean > 0),
                            cores = 4)

Sn_ln_fS_region <- update(Sn_lognorm_fragSize, 
                          formula. = ~ c.lfs * continent8 + (c.lfs | dataset_label),
                          newdata = frag %>% filter(S_n_mean > 0),
                          cores = 4)

Scov_ln_fS_region <- update(Scov_lognorm_fragSize, 
                            formula. = ~ c.lfs * continent8 + (c.lfs | dataset_label),
                            newdata = frag %>% filter(S_cov_mean > 0),
                            cores = 4)

Schao_ln_fS_region <- update(S_chao_lognorm_fragSize, 
                             formula. = ~ c.lfs * continent8 + (c.lfs | dataset_label),
                             newdata = frag %>% filter(S_chao_mean > 0),
                             cores = 4)

S_PIE_ln_fS_region <- update(S_PIE_lognorm_fragSize, 
                             formula. = ~ c.lfs * continent8 + (c.lfs | dataset_label),
                             newdata = frag %>% filter(S_PIE_mean > 0),
                             cores = 4)

N_std_ln_fS_region <- update(Nstd_lognorm_fragSize, 
                             formula. = ~ c.lfs * continent8 + (c.lfs | dataset_label),
                             newdata = frag,
                             cores = 4)  

# compare the model fits (with versus without interactions
Sstd_lognorm_fragSize <- add_criterion(Sstd_lognorm_fragSize, criterion =  'loo') 
Sstd_ln_fS_matrix <- add_criterion(Sstd_ln_fS_matrix, criterion = 'loo')
Sstd_ln_fS_taxa <- add_criterion(Sstd_ln_fS_taxa, criterion = 'loo')
Sstd_ln_fS_tsf <- add_criterion(Sstd_ln_fS_tsf, criterion = 'loo')
Sstd_ln_fS_biome <- add_criterion(Sstd_ln_fS_biome, criterion = 'loo')
Sstd_ln_fS_region <- add_criterion(Sstd_ln_fS_region, criterion = 'loo')

loo_compare(Sstd_lognorm_fragSize,
            Sstd_ln_fS_matrix,
            Sstd_ln_fS_taxa,
            Sstd_ln_fS_tsf,
            Sstd_ln_fS_biome,
            Sstd_ln_fS_region)

S_PIE_lognorm_fragSize <- add_criterion(S_PIE_lognorm_fragSize, criterion =  'loo') 
S_PIE_ln_fS_matrix <- add_criterion(S_PIE_ln_fS_matrix, criterion = 'loo')
S_PIE_ln_fS_taxa <- add_criterion(S_PIE_ln_fS_taxa, criterion = 'loo')
S_PIE_ln_fS_tsf <- add_criterion(S_PIE_ln_fS_tsf, criterion = 'loo')
S_PIE_ln_fS_biome <- add_criterion(S_PIE_ln_fS_biome, criterion = 'loo')
S_PIE_ln_fS_region <- add_criterion(S_PIE_ln_fS_region, criterion = 'loo')

loo_compare(S_PIE_lognorm_fragSize,
            S_PIE_ln_fS_matrix,
            S_PIE_ln_fS_taxa,
            S_PIE_ln_fS_tsf,
            S_PIE_ln_fS_biome,
            S_PIE_ln_fS_region)

Nstd_lognorm_fragSize <- add_criterion(Nstd_lognorm_fragSize, criterion =  'loo') 
N_std_ln_fS_matrix <- add_criterion(N_std_ln_fS_matrix, criterion = 'loo')
N_std_ln_fS_taxa <- add_criterion(N_std_ln_fS_taxa, criterion = 'loo')
N_std_ln_fS_tsf <- add_criterion(N_std_ln_fS_tsf, criterion = 'loo')
N_std_ln_fS_biome <- add_criterion(N_std_ln_fS_biome, criterion = 'loo')
N_std_ln_fS_region <- add_criterion(N_std_ln_fS_region, criterion = 'loo')

loo_compare(Nstd_lognorm_fragSize,
            N_std_ln_fS_matrix,
            N_std_ln_fS_taxa,
            N_std_ln_fS_tsf,
            N_std_ln_fS_biome,
            N_std_ln_fS_region)


save(Sstd_ln_fS_matrix, Sstd_ln_fS_taxa, Sstd_ln_fS_tsf, Sstd_ln_fS_biome, Sstd_ln_fS_region,
     Sn_ln_fS_matrix, Sn_ln_fS_taxa, Sn_ln_fS_tsf, Sn_ln_fS_biome,
     Scov_ln_fS_matrix, Scov_ln_fS_taxa, Scov_ln_fS_tsf, Scov_ln_fS_biome,
     Schao_ln_fS_matrix, Schao_ln_fS_taxa, Schao_ln_fS_tsf, Schao_ln_fS_biome,
     S_PIE_ln_fS_matrix, S_PIE_ln_fS_taxa, S_PIE_ln_fS_tsf, S_PIE_ln_fS_biome, S_PIE_ln_fS_region,
     N_std_ln_fS_matrix, N_std_ln_fS_taxa, N_std_ln_fS_tsf, N_std_ln_fS_biome, N_std_ln_fS_region,
     Sstd_lognorm_fragSize, S_PIE_lognorm_fragSize, Nstd_lognorm_fragSize,
     file = '~/Dropbox/1current/fragmentation_synthesis/results/fragSize_interactions_ref_revision.Rdata')


