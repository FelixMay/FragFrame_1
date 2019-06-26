# get study-level posterior samples from model with only study-level slope variation: beta diversity
library(tidyverse)
library(brms)
library(ggridges)

load('~/Dropbox/1current/fragmentation_synthesis/results/jtu_z1i_fS_noPairGroup-5282398.Rdata')
load('~/Dropbox/1current/fragmentation_synthesis/results/rtu_z1i_fS_noPairGroup-5282404.Rdata')

load('~/Dropbox/1current/fragmentation_synthesis/results/Jne_zi_fragSize_noPairGroup.Rdata')
load('~/Dropbox/1current/fragmentation_synthesis/results/Rne_zi_fragSize_noPairGroup.Rdata')

frag_beta <- read_csv('~/Dropbox/Frag Database (new)/files_datapaper/Analysis/2_betapart_frag_fcont_10_mabund_as_is.csv')

# get the metadata
meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

# check names
meta_labels <- meta %>% distinct(dataset_label)

meta_labels %>% 
  filter(!dataset_label %in% Jtu_z1i_fS$data$dataset_label) %>% 
  distinct(dataset_label)

Jtu_z1i_fS$data %>% 
  filter(!dataset_label %in% meta_labels$dataset_label) %>% 
  distinct(dataset_label)

# change metadata labels (as the ones in frag were used for the model fitting)
meta <- meta %>% 
  mutate(dataset_label = as.character(dataset_label),
         dataset_label = ifelse(dataset_label=='Brosi_2009', 'Brosi_2007', dataset_label),
         dataset_label = ifelse(dataset_label=='delaSancha_2014', 'DeLaSancha_2014', dataset_label),
         dataset_label = ifelse(dataset_label=='deSouza_1994', 'DeSouza_1994', dataset_label),
         dataset_label = ifelse(dataset_label=='Dominguez-Haydar_2011', 'Dominquez-Haydar_2011', dataset_label),
         dataset_label = ifelse(dataset_label=='Guadagnin_2005', 'Gaudagnin_2005', dataset_label),
         dataset_label = ifelse(dataset_label=='Raheem_2009', 'Raheen_2009', dataset_label),
         dataset_label = ifelse(dataset_label=='Silveira_2015', 'Silviera_2015', dataset_label),
         dataset_label = ifelse(dataset_label=='Telleria_1995', 'Tellbera_1995', dataset_label),
         dataset_label = ifelse(dataset_label=='Vulinec_2008', 'Vulineci_2008', dataset_label),
         dataset_label = ifelse(dataset_label=='Sridhar_2008', 'Sridihar_2008', dataset_label))


frag_beta <- frag_beta %>% 
  # centre covariate before fitting
  mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))

frag_beta <- left_join(frag_beta, 
                       meta, by = 'dataset_label')


# study-levels (no studies are missing N_std)
study_levels <- Jtu_z1i_fS$data %>% 
  as_tibble() %>% 
  distinct(dataset_label) %>% 
  mutate(level = dataset_label) %>%
  nest(level) 

study_sample_posterior <- study_levels %>%
  mutate(Jtu = purrr::map(data, ~posterior_samples(Jtu_z1i_fS, 
                                                     pars = paste('r_dataset_label[', as.character(.x$level), ',cl10ra]', sep=''),
                                                     exact = TRUE,
                                                     subset = floor(runif(n = 1000,
                                                                          min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         Rtu = purrr::map(data, ~posterior_samples(Rtu_z1i_fS, 
                                                    pars = paste('r_dataset_label[', as.character(.x$level), ',cl10ra]', sep=''),
                                                    exact = TRUE,
                                                    subset = floor(runif(n = 1000,
                                                                         min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         Jne = purrr::map(data, ~posterior_samples(Jne_zi_fragSize,
                                                  pars = paste('r_dataset_label[', as.character(.x$level), ',cl10ra]', sep=''),
                                                  exact = TRUE,
                                                  subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()),
         Rne = purrr::map(data, ~posterior_samples(Rne_zi_fragSize,
                                                     pars = paste('r_dataset_label[', as.character(.x$level), ',cl10ra]', sep=''),
                                                     exact = TRUE,
                                                     subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()))


Jtu_fixef <- fixef(Jtu_z1i_fS)
Rtu_fixef <- fixef(Rtu_z1i_fS)
Jne_fixef <- fixef(Jne_zi_fragSize)
Rne_fixef <- fixef(Rne_zi_fragSize)

Jtu_posterior <- study_sample_posterior  %>% 
  select(-data) %>% 
  unnest(Jtu) %>% 
  mutate(response = 'Jtu',
         Jtu_global_slope = Jtu_fixef['cl10ra','Estimate'],
         Jtu_upper_slope = Jtu_fixef['cl10ra','Q97.5'],
         Jtu_lower_slope = Jtu_fixef['cl10ra','Q2.5']) %>% 
  left_join(meta, 
            by = 'dataset_label')

Jtu_posterior$time.since.fragmentation <- factor(Jtu_posterior$time.since.fragmentation,
                                                  levels = c('Recent (less than 20 years)',
                                                             'Intermediate (20-100 years)',
                                                             'long (100+ years)'),
                                                  labels = c('< 20 years',
                                                             '20-100 years',
                                                             '> 100 years'))

Jtu_posterior$Matrix.category <- factor(Jtu_posterior$Matrix.category,
                                         levels = c('light filter', 'intermediate', 'harsh filter'),
                                         labels = c('Light', 'Intermediate', 'Harsh'))

Jtu_posterior$biome <- factor(Jtu_posterior$biome,
                               levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                               labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))
Jtu_posterior$taxa <- factor(Jtu_posterior$taxa,
                              levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                              labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

Jtu_posterior <- Jtu_posterior %>% 
  unite(col = 'frag_matrix', c(sphere.fragment, sphere.matrix),
        remove = F, sep = ', ')

# rtu
Rtu_posterior <- study_sample_posterior  %>% 
  select(-data) %>% 
  unnest(Rtu) %>% 
  mutate(response = 'Rtu',
         Rtu_global_slope = Rtu_fixef['cl10ra','Estimate'],
         Rtu_upper_slope = Rtu_fixef['cl10ra','Q97.5'],
         Rtu_lower_slope = Rtu_fixef['cl10ra','Q2.5']) %>% 
  left_join(meta, 
            by = 'dataset_label')

Rtu_posterior$time.since.fragmentation <- factor(Rtu_posterior$time.since.fragmentation,
                                                 levels = c('Recent (less than 20 years)',
                                                            'Intermediate (20-100 years)',
                                                            'long (100+ years)'),
                                                 labels = c('< 20 years',
                                                            '20-100 years',
                                                            '> 100 years'))

Rtu_posterior$Matrix.category <- factor(Rtu_posterior$Matrix.category,
                                        levels = c('light filter', 'intermediate', 'harsh filter'),
                                        labels = c('Light', 'Intermediate', 'Harsh'))

Rtu_posterior$biome <- factor(Rtu_posterior$biome,
                              levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                              labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))
Rtu_posterior$taxa <- factor(Rtu_posterior$taxa,
                             levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                             labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

Rtu_posterior <- Rtu_posterior %>% 
  unite(col = 'frag_matrix', c(sphere.fragment, sphere.matrix),
        remove = F, sep = ', ')

# Jne
Jne_posterior <- study_sample_posterior  %>% 
  select(-data) %>% 
  unnest(Jne) %>% 
  mutate(response = 'Jne',
         Jne_global_slope = Jne_fixef['cl10ra','Estimate'],
         Jne_upper_slope = Jne_fixef['cl10ra','Q97.5'],
         Jne_lower_slope = Jne_fixef['cl10ra','Q2.5']) %>% 
  left_join(meta, 
            by = 'dataset_label')

Jne_posterior$time.since.fragmentation <- factor(Jne_posterior$time.since.fragmentation,
                                                 levels = c('Recent (less than 20 years)',
                                                            'Intermediate (20-100 years)',
                                                            'long (100+ years)'),
                                                 labels = c('< 20 years',
                                                            '20-100 years',
                                                            '> 100 years'))

Jne_posterior$Matrix.category <- factor(Jne_posterior$Matrix.category,
                                        levels = c('light filter', 'intermediate', 'harsh filter'),
                                        labels = c('Light', 'Intermediate', 'Harsh'))

Jne_posterior$biome <- factor(Jne_posterior$biome,
                              levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                              labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))
Jne_posterior$taxa <- factor(Jne_posterior$taxa,
                             levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                             labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

Jne_posterior <- Jne_posterior %>% 
  unite(col = 'frag_matrix', c(sphere.fragment, sphere.matrix),
        remove = F, sep = ', ')

# Rne
Rne_posterior <- study_sample_posterior  %>% 
  select(-data) %>% 
  unnest(Rne) %>% 
  mutate(response = 'Rne',
         Rne_global_slope = Rne_fixef['cl10ra','Estimate'],
         Rne_upper_slope = Rne_fixef['cl10ra','Q97.5'],
         Rne_lower_slope = Rne_fixef['cl10ra','Q2.5']) %>% 
  left_join(meta, 
            by = 'dataset_label')

Rne_posterior$time.since.fragmentation <- factor(Rne_posterior$time.since.fragmentation,
                                                 levels = c('Recent (less than 20 years)',
                                                            'Intermediate (20-100 years)',
                                                            'long (100+ years)'),
                                                 labels = c('< 20 years',
                                                            '20-100 years',
                                                            '> 100 years'))

Rne_posterior$Matrix.category <- factor(Rne_posterior$Matrix.category,
                                        levels = c('light filter', 'intermediate', 'harsh filter'),
                                        labels = c('Light', 'Intermediate', 'Harsh'))

Rne_posterior$biome <- factor(Rne_posterior$biome,
                              levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                              labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))
Rne_posterior$taxa <- factor(Rne_posterior$taxa,
                             levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                             labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

Rne_posterior <- Rne_posterior %>% 
  unite(col = 'frag_matrix', c(sphere.fragment, sphere.matrix),
        remove = F, sep = ', ')

