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

# # set the reference levels for the categorical covariates of interest
frag %>% distinct(Matrix.category)
frag$Matrix.category <- factor(frag$Matrix.category, levels = c('light filter', 'intermediate', 'harsh filter'))
frag %>% distinct(time.since.fragmentation)
frag$time.since.fragmentation <- factor(frag$time.since.fragmentation,
                                         levels = c("Recent (less than 20 years)", "Intermediate (20-100 years)", "long (100+ years)"))


# two-way interactions: matrix permeability first
Sstd2_ln_fS_matrix <- update(Sstd2_lognorm_fragSize,
                             formula. = ~ c.lfs * Matrix.category + (c.lfs | dataset_label),
                             newdata = frag, cores = 4)

Sn_ln_fS_matrix <- update(Sn_lognorm_fragSize, 
                          formula. = ~ c.lfs * Matrix.category + (c.lfs | dataset_label),
                          newdata = frag, cores = 4)

Scov_ln_fS_matrix <- update(Scov_lognorm_fragSize, 
                          formula. = ~ c.lfs * Matrix.category + (c.lfs | dataset_label),
                          newdata = frag, cores = 4)

Schao_ln_fS_matrix <- update(S_chao_lognorm_fragSize, 
                             formula. = ~ c.lfs * Matrix.category + (c.lfs | dataset_label),
                             newdata = frag, cores = 4)

S_PIE_ln_fS_matrix <- update(S_PIE_lognorm_fragSize, 
                             formula. = ~ c.lfs * Matrix.category + (c.lfs | dataset_label),
                             newdata = frag, cores = 4)

N_std_ln_fS_matrix <- update(Nstd_lognorm_fragSize, 
                             formula. = ~ c.lfs * Matrix.category + (c.lfs | dataset_label),
                             newdata = frag, cores = 4)  

# repeat for taxa
Sstd2_ln_fS_taxa <- update(Sstd2_lognorm_fragSize, 
                          formula. = ~ c.lfs * taxa + (c.lfs | dataset_label),
                          newdata = frag, cores = 4)

Sn_ln_fS_taxa <- update(Sn_lognorm_fragSize, 
                          formula. = ~ c.lfs * taxa + (c.lfs | dataset_label),
                          newdata = frag, cores = 4)

Scov_ln_fS_taxa <- update(Scov_lognorm_fragSize, 
                            formula. = ~ c.lfs * taxa + (c.lfs | dataset_label),
                            newdata = frag, cores = 4)

Schao_ln_fS_taxa <- update(S_chao_lognorm_fragSize, 
                             formula. = ~ c.lfs * taxa + (c.lfs | dataset_label),
                             newdata = frag, cores = 4)

S_PIE_ln_fS_taxa <- update(S_PIE_lognorm_fragSize, 
                             formula. = ~ c.lfs * taxa + (c.lfs | dataset_label),
                             newdata = frag, cores = 4)

N_std_ln_fS_taxa <- update(Nstd_lognorm_fragSize, 
                             formula. = ~ c.lfs * taxa + (c.lfs | dataset_label),
                             newdata = frag, cores = 4)  

# repeat for time since fragmentation
Sstd2_ln_fS_tsf <- update(Sstd2_lognorm_fragSize, 
                           formula. = ~ c.lfs * time.since.fragmentation + (c.lfs | dataset_label),
                           newdata = frag, cores = 4)

Sn_ln_fS_tsf <- update(Sn_lognorm_fragSize, 
                        formula. = ~ c.lfs * time.since.fragmentation + (c.lfs | dataset_label),
                        newdata = frag, cores = 4)

Scov_ln_fS_tsf <- update(Scov_lognorm_fragSize, 
                          formula. = ~ c.lfs * time.since.fragmentation + (c.lfs | dataset_label),
                          newdata = frag, cores = 4)

Schao_ln_fS_tsf <- update(S_chao_lognorm_fragSize, 
                           formula. = ~ c.lfs * time.since.fragmentation + (c.lfs | dataset_label),
                           newdata = frag, cores = 4)

S_PIE_ln_fS_tsf <- update(S_PIE_lognorm_fragSize, 
                           formula. = ~ c.lfs * time.since.fragmentation + (c.lfs | dataset_label),
                           newdata = frag, cores = 4)

N_std_ln_fS_tsf <- update(Nstd_lognorm_fragSize, 
                           formula. = ~ c.lfs * time.since.fragmentation + (c.lfs | dataset_label),
                           newdata = frag, cores = 4)  

# repeat for biome (remove the single wetland study)
Sstd2_ln_fS_biome <- update(Sstd2_lognorm_fragSize, 
                          formula. = ~ c.lfs * biome + (c.lfs | dataset_label),
                          newdata = frag, #%>% filter(biome!='wetland'),
                          cores = 4)

Sn_ln_fS_biome <- update(Sn_lognorm_fragSize, 
                       formula. = ~ c.lfs * biome + (c.lfs | dataset_label),
                       newdata = frag, #%>% filter(biome!='wetland'), 
                       cores = 4)

Scov_ln_fS_biome <- update(Scov_lognorm_fragSize, 
                         formula. = ~ c.lfs * biome + (c.lfs | dataset_label),
                         newdata = frag,# %>% filter(biome!='wetland'), 
                         cores = 4)

Schao_ln_fS_biome <- update(S_chao_lognorm_fragSize, 
                          formula. = ~ c.lfs * biome + (c.lfs | dataset_label),
                          newdata = frag,# %>% filter(biome!='wetland'), 
                          cores = 4)

S_PIE_ln_fS_biome <- update(S_PIE_lognorm_fragSize, 
                          formula. = ~ c.lfs * biome + (c.lfs | dataset_label),
                          newdata = frag, #%>% filter(biome!='wetland'), 
                          cores = 4)

N_std_ln_fS_biome <- update(Nstd_lognorm_fragSize, 
                            formula. = ~ c.lfs * biome + (c.lfs | dataset_label),
                            newdata = frag,# %>% filter(biome!='wetland'), 
                            cores = 4)

# I initially fit the biome models without the wetland study
save(Sstd2_ln_fS_matrix, Sstd2_ln_fS_taxa, Sstd2_ln_fS_tsf, Sstd2_ln_fS_biome,
     Sn_ln_fS_matrix, Sn_ln_fS_taxa, Sn_ln_fS_tsf, Sn_ln_fS_biome,
     Scov_ln_fS_matrix, Scov_ln_fS_taxa, Scov_ln_fS_tsf, Scov_ln_fS_biome,
     Schao_ln_fS_matrix, Schao_ln_fS_taxa, Schao_ln_fS_tsf, Schao_ln_fS_biome,
     S_PIE_ln_fS_matrix, S_PIE_ln_fS_taxa, S_PIE_ln_fS_tsf, S_PIE_ln_fS_biome,
     N_std_ln_fS_matrix, N_std_ln_fS_taxa, N_std_ln_fS_tsf, N_std_ln_fS_biome,
     file = '~/Dropbox/1current/fragmentation_synthesis/results/fragSize_interactions.Rdata')

# refit biome models with wetland
save(Sstd2_ln_fS_biome,
     Sn_ln_fS_biome,
     Scov_ln_fS_biome,
     Schao_ln_fS_biome,
     S_PIE_ln_fS_biome,
     N_std_ln_fS_biome,
     file = '~/Dropbox/1current/fragmentation_synthesis/results/fragSize_biome_wetland.Rdata')

# compare the model fits (with versus without interactions)
waic(Sstd2_lognorm_fragSize, 
     Sstd2_ln_fS_matrix,
     Sstd2_ln_fS_taxa,
     Sstd2_ln_fS_tsf,
     Sstd2_ln_fS_biome)

waic(Sn_lognorm_fragSize,
     Sn_ln_fS_matrix,
     Sn_ln_fS_taxa,
     Sn_ln_fS_tsf,
     Sn_ln_fS_biome)

waic(Scov_lognorm_fragSize,
     Scov_ln_fS_matrix,
     Scov_ln_fS_taxa,
     Scov_ln_fS_tsf,
     Scov_ln_fS_biome)

waic(S_chao_lognorm_fragSize,
     Schao_ln_fS_matrix,
     Schao_ln_fS_taxa,
     Schao_ln_fS_tsf,
     Schao_ln_fS_biome)

waic(S_PIE_lognorm_fragSize,
     S_PIE_ln_fS_matrix,
     S_PIE_ln_fS_taxa,
     S_PIE_ln_fS_tsf,
     S_PIE_ln_fS_biome)

waic(Nstd_lognorm_fragSize,
     N_std_ln_fS_matrix,
     N_std_ln_fS_taxa,
     N_std_ln_fS_tsf,
     N_std_ln_fS_biome)
