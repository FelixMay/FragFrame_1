# execute 0_init_dirs_load_packages.R for packages and directory

# get posterior samples from model for overall trend and study-level slopes: this time the z-score regressions

# load model fits and the data
load(paste0(path2wd, 'intermediate_results/fragSize_z_score_ref.Rdata'))

frag <- read_csv(paste0(path2wd, 'intermediate_results/2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

# load the meta data
meta <- read.csv(paste0(path2wd, 'data/new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')

# get posterior distribution of global estimates
frag_global <- tibble(
  Sstd_global = posterior_samples(z_Sstd_studT_fragSize,
                                  pars = 'b_c.lfs',
                                  exact_match = TRUE,
                                  subset = floor(runif(n = 1000,
                                                       min = 1, max = 2000))) %>% unlist() %>% as.numeric(),
  Scov_global = posterior_samples(z_Scov_studT_fragSize,
                                  pars = 'b_c.lfs',
                                  exact_match = TRUE,
                                  subset = floor(runif(n = 1000,
                                                       min = 1, max = 2000))) %>% unlist() %>% as.numeric(),
  Sn_global = posterior_samples(z_Sn_studT_fragSize,
                                pars = 'b_c.lfs',
                                exact_match = TRUE,
                                subset = floor(runif(n = 1000,
                                                     min = 1, max = 2000))) %>% unlist() %>% as.numeric(),
  Schao_global = posterior_samples(z_S_chao_studT_fragSize,
                                   pars = 'b_c.lfs',
                                   exact_match = TRUE,
                                   subset = floor(runif(n = 1000,
                                                        min = 1, max = 2000))) %>% unlist() %>% as.numeric(),
  S_PIE_global = posterior_samples(z_S_PIE_studT_fragSize,
                                   pars = 'b_c.lfs',
                                   exact_match = TRUE,
                                   subset = floor(runif(n = 1000,
                                                        min = 1, max = 2000))) %>% unlist() %>% as.numeric()
)

# study-levels 
study_levels <- z_Sstd_studT_fragSize$data %>% 
  as_tibble() %>% 
  distinct(dataset_label) %>% 
  mutate(level = dataset_label) %>%
  nest(level) 

study_sample_posterior <- study_levels %>%
  mutate(S_std = purrr::map(data, ~posterior_samples(z_Sstd_studT_fragSize, 
                                                     pars = paste('r_dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                     exact = TRUE,
                                                     subset = floor(runif(n = 1000,
                                                                          min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         Scov = purrr::map(data, ~posterior_samples(z_Scov_studT_fragSize, 
                                                    pars = paste('r_dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                    exact = TRUE,
                                                    subset = floor(runif(n = 1000,
                                                                         min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         Sn = purrr::map(data, ~posterior_samples(z_Sn_studT_fragSize,
                                                  pars = paste('r_dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                  exact = TRUE,
                                                  subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()),
         Schao = purrr::map(data, ~posterior_samples(z_S_chao_studT_fragSize,
                                                     pars = paste('r_dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                     exact = TRUE,
                                                     subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()),
         S_PIE = purrr::map(data, ~posterior_samples(z_S_PIE_studT_fragSize,
                                                     pars = paste('r_dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                     exact = TRUE,
                                                     subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric())
  )


Sstd_posterior <- study_sample_posterior  %>% 
  select(-data) %>% 
  unnest(S_std) %>% 
  mutate(response = 'z_Sstd',
         Sstd_global = rep(frag_global$Sstd_global, times = n_distinct(dataset_label))) %>% 
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
  mutate(response = 'z_S_PIE',
         S_PIE_global = rep(frag_global$S_PIE_global, times = n_distinct(dataset_label))) %>% 
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

Spie_posterior$biome <- factor(Spie_posterior$biome,
                               levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                               labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))
Spie_posterior$taxa <- factor(Spie_posterior$taxa,
                              levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                              labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

