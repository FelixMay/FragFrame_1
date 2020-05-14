# need to execute 0_init_dirs_load_packages.R first

# code to fit models with 2-way interaction for fragmentation synthesis

frag <- read_csv(paste0(path2wd, 'intermediate_results/2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

meta <- read_delim(paste0(path2wd, 'data/new_meta_2_merge.csv'),  delim =';') %>% 
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
load(paste0(path2wd, 'main_results/fragSize_ref.Rdata'))


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


# save
save(Sstd_ln_fS_matrix, Sstd_ln_fS_taxa, Sstd_ln_fS_tsf, Sstd_ln_fS_region,
     file = paste0(path2wd, 'extended_data_figs_tabs/Sstd_interactionModels.Rdata'))

save(N_std_ln_fS_matrix, N_std_ln_fS_taxa, N_std_ln_fS_tsf, N_std_ln_fS_region,
     file = paste0(path2wd, 'extended_data_figs_tabs/Nstd_interactionModels.Rdata'))

save(S_PIE_ln_fS_matrix, S_PIE_ln_fS_region, S_PIE_ln_fS_tsf, S_PIE_ln_fS_taxa,
     file = paste0(path2wd, 'extended_data_figs_tabs/S_PIE_interactionModels.Rdata'))
