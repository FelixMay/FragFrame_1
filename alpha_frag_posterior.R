# get study-level posterior samples
library(tidyverse)
library(brms)
library(ggridges)

load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms.Rdata')

# get the metadata (I want to group the posteriors)
meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

# higher-levels (no studies are missing N_std)
high_level <- Nstd_lognorm_fragSize4$data %>% 
  as_tibble() %>% 
  distinct(cbmtt) %>% 
  mutate(level = cbmtt) %>%
  nest(level)

# study-levels (no studies are missing N_std)
study_levels <- Nstd_lognorm_fragSize$data %>% 
  as_tibble() %>% 
  distinct(dataset_label) %>% 
  mutate(level = dataset_label) %>%
  nest(level)

cbmtt_posterior <- high_level %>%
  mutate(S_std = purrr::map(data, possibly(~posterior_samples(Sstd2_lognorm_fragSize4,
                                                     pars = paste0('r_cbmtt[', as.character(.x$level), ',c.lfs]'),
                                                     exact_match = TRUE,
                                                     subset = floor(runif(n = 1000,
                                                                          min = 1, max = 2000))) %>% unlist() %>% as.numeric(),
                                           otherwise = NULL)),
         Scov = purrr::map(data, ~posterior_samples(Scov_lognorm_fragSize4,
                                                    pars = paste('r_cbmtt[', as.character(.x$level), ',c.lfs]', sep=''),
                                                    exact_match = TRUE,
                                                    subset = floor(runif(n = 1000,
                                                                         min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         Sn = purrr::map(data, ~posterior_samples(Sn_lognorm_fragSize4,
                                                  pars = paste('r_cbmtt[', as.character(.x$level), ',c.lfs]', sep=''),
                                                  exact_match = TRUE,
                                                  subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()),
         Schao = purrr::map(data, ~posterior_samples(S_chao_lognorm_fragSize4,
                                                     pars = paste('r_cbmtt[', as.character(.x$level), ',c.lfs]', sep=''),
                                                     exact_match = TRUE,
                                                     subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()),
         S_PIE = purrr::map(data, ~posterior_samples(S_PIE_lognorm_fragSize4,
                                                     pars = paste('r_cbmtt[', as.character(.x$level), ',c.lfs]', sep=''),
                                                     exact_match = TRUE,
                                                     subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()),
         Nstd = purrr::map(data, ~posterior_samples(Nstd_lognorm_fragSize4,
                                                    pars = paste('r_cbmtt[', as.character(.x$level), ',c.lfs]', sep=''),
                                                    exact_match = TRUE,
                                                    subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()))

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
                                                             '20-100 years)',
                                                             '> 100 years)'))

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
ggplot() +
  facet_grid(continent ~ ., scale = 'free') +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(xmin = Sstd_lower_slope, xmax = Sstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = Sstd_posterior,
                      aes(x = S_std + unique(Sstd_global_slope), 
                          y = time.since.fragmentation,
                          fill = taxa
                      ),
                      scale = 1, alpha = 0.2,
                      linetype = 0) +
  # scale_fill_viridis_d(name = 'Taxa') +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Time since fragmentation',
       x = expression(paste('Study-level ', S[std], ' change'))#,
       # tag = 'A'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.85, 0.8)) #+

ggplot() +
  facet_grid(~ taxa, scales = 'free') +
  geom_rect(data = Spie_posterior %>% distinct(S_PIE_lower_slope, S_PIE_upper_slope),
            aes(xmin = S_PIE_lower_slope, xmax = S_PIE_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = Spie_posterior %>% filter(biome=='forest'),
                      aes(x = S_PIE + unique(S_PIE_global_slope), 
                          y = time.since.fragmentation,#interaction(Matrix.category, time.since.fragmentation),
                          fill = Matrix.category
                      ),
                      scale = 0.9, alpha = 0.6,
                      linetype = 0) +
  # scale_linetype_manual(name = 'Realm', values = c('Marine' = 0, 'Terrestrial' = 1, 'Freshwater' = 3)) +
  # scale_fill_manual(name = 'Taxa', values = taxa_col) +
  geom_vline(data = Spie_posterior,
             aes(xintercept = S_PIE_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Time since fragmentation & matrix category',
       x = expression(paste('Study-level ', S[PIE], ' change'))#,
       # tag = 'A'
  ) +
  # scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.7, 0.2),
        # legend.direction = 'horizontal',
        legend.background = element_blank()) #+

