# get study-level posterior samples from model with only study-level slope variation
library(tidyverse)
library(brms)
library(ggridges)

load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms.Rdata')

# get the metadata (I want to group the posteriors)
meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)


# study-levels (no studies are missing N_std)
study_levels <- Nstd_lognorm_fragSize$data %>% 
  as_tibble() %>% 
  distinct(dataset_label) %>% 
  mutate(level = dataset_label) %>%
  nest(level) 

study_sample_posterior <- study_levels %>%
  mutate(S_std = purrr::map(data, ~posterior_samples(Sstd2_lognorm_fragSize, 
                                                     pars = paste('r_dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                     exact = TRUE,
                                                     subset = floor(runif(n = 1000,
                                                                          min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         Scov = purrr::map(data, ~posterior_samples(Scov_lognorm_fragSize, 
                                                    pars = paste('r_dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                    exact = TRUE,
                                                    subset = floor(runif(n = 1000,
                                                                         min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         Sn = purrr::map(data, ~posterior_samples(Sn_lognorm_fragSize,
                                                  pars = paste('r_dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                  exact = TRUE,
                                                  subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()),
         Schao = purrr::map(data, ~posterior_samples(S_chao_lognorm_fragSize,
                                                  pars = paste('r_dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                  exact = TRUE,
                                                  subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()),
         S_PIE = purrr::map(data, ~posterior_samples(S_PIE_lognorm_fragSize,
                                                     pars = paste('r_dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                     exact = TRUE,
                                                     subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()),
         Nstd = purrr::map(data, ~posterior_samples(Nstd_lognorm_fragSize,
                                                    pars = paste('r_dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                    exact = TRUE,
                                                    subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()))


Sstd2_lognorm_fragSize_fixef <- fixef(Sstd2_lognorm_fragSize)
S_PIE_lognorm_fragSize_fixef <- fixef(S_PIE_lognorm_fragSize)
Scov_lognorm_fragSize_fixef <- fixef(Scov_lognorm_fragSize)
Sn_lognorm_fragSize_fixef <- fixef(Sn_lognorm_fragSize)
Schao_lognorm_fragSize_fixef <- fixef(S_chao_lognorm_fragSize)
Nstd_lognorm_fragSize_fixef <- fixef(Nstd_lognorm_fragSize)

Sstd_posterior <- study_sample_posterior  %>% 
  select(-data) %>% 
  unnest(S_std) %>% 
  mutate(response = 'Sstd',
         Sstd_global_slope = Sstd2_lognorm_fragSize_fixef['c.lfs','Estimate'],
         Sstd_upper_slope = Sstd2_lognorm_fragSize_fixef['c.lfs','Q97.5'],
         Sstd_lower_slope = Sstd2_lognorm_fragSize_fixef['c.lfs','Q2.5']) %>% 
  left_join(meta, 
            by = 'dataset_label')

Sstd_posterior$time.since.fragmentation <- factor(Sstd_posterior$time.since.fragmentation,
                                                  levels = c('Recent (less than 20 years)',
                                                             'Intermediate (20-100 years)',
                                                             'long (100+ years)'),
                                                  labels = c('< 20 years',
                                                             '20-100 years',
                                                             '> 100 years'))

Sstd_posterior$Matrix.category <- factor(Sstd_posterior$Matrix.category,
                                                  levels = c('light filter', 'intermediate', 'harsh filter'),
                                                  labels = c('Light',
                                                             'Intermediate',
                                                             'Harsh'))

Sstd_posterior$biome <- factor(Sstd_posterior$biome,
                                         levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                                         labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))
Sstd_posterior$taxa <- factor(Sstd_posterior$taxa,
                               levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                               labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

Sstd_posterior <- Sstd_posterior %>% 
  unite(col = 'frag_matrix', c(sphere.fragment, sphere.matrix),
        remove = F, sep = ', ')

Spie_posterior <- study_sample_posterior  %>% 
  select(-data) %>% 
  unnest(S_PIE) %>% 
  mutate(response = 'S_PIE',
         S_PIE_global_slope = S_PIE_lognorm_fragSize_fixef['c.lfs','Estimate'],
         S_PIE_upper_slope = S_PIE_lognorm_fragSize_fixef['c.lfs','Q97.5'],
         S_PIE_lower_slope = S_PIE_lognorm_fragSize_fixef['c.lfs','Q2.5']) %>% 
  left_join(meta, 
            by = 'dataset_label')

Spie_posterior$time.since.fragmentation <- factor(Spie_posterior$time.since.fragmentation,
                                                  levels = c('Recent (less than 20 years)',
                                                             'Intermediate (20-100 years)',
                                                             'long (100+ years)'),
                                                  labels = c('< 20 years',
                                                             '20-100 years)',
                                                             '> 100 years)'))
Spie_posterior$Matrix.category <- factor(Spie_posterior$Matrix.category,
                                         levels = c('light filter', 'intermediate', 'harsh filter'),
                                         labels = c('Light',
                                                    'Intermediate',
                                                    'Harsh'))

