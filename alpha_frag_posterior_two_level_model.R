# get  posterior samples from model with two levels: cbmtt AND study-level slope variation
library(tidyverse)
library(brms)
library(ggridges)

load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms.Rdata')

# get the metadata (I want to group the posteriors)
meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

# cbmtt levels
high_level <- Nstd_lognorm_fragSize4$data %>%
  as_tibble() %>%
  distinct(cbmtt) %>%
  mutate(level = cbmtt,
         level = gsub(" ", ".", level)) %>%
  nest(level)

# study-levels (no studies are missing N_std)
study_levels <- Nstd_lognorm_fragSize4$data %>% 
  as_tibble() %>% 
  distinct(cbmtt, dataset_label) %>% 
  unite(level, c(cbmtt,dataset_label), remove = F) %>%
  mutate(level = gsub(" ", ".", level)) %>% 
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
  mutate(S_std = purrr::map(data, ~posterior_samples(Sstd2_lognorm_fragSize4, 
                                                     pars = paste('r_cbmtt:dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                     exact = TRUE,
                                                     subset = floor(runif(n = 1000,
                                                                          min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         Scov = purrr::map(data, ~posterior_samples(Scov_lognorm_fragSize4, 
                                                    pars = paste('r_cbmtt:dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                    exact = TRUE,
                                                    subset = floor(runif(n = 1000,
                                                                         min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         Sn = purrr::map(data, ~posterior_samples(Sn_lognorm_fragSize4,
                                                  pars = paste('r_cbmtt:dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                  exact = TRUE,
                                                  subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()),
         Schao = purrr::map(data, ~posterior_samples(S_chao_lognorm_fragSize4,
                                                     pars = paste('r_cbmtt:dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                     exact = TRUE,
                                                     subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()),
         S_PIE = purrr::map(data, ~posterior_samples(S_PIE_lognorm_fragSize4,
                                                     pars = paste('r_cbmtt:dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                     exact = TRUE,
                                                     subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()),
         Nstd = purrr::map(data, ~posterior_samples(Nstd_lognorm_fragSize4,
                                                    pars = paste('r_cbmtt:dataset_label[', as.character(.x$level), ',c.lfs]', sep=''),
                                                    exact = TRUE,
                                                    subset = floor(runif(n = 1000, 1, max = 2000))) %>%  unlist() %>%  as.numeric()))


Sstd2_lognorm_fragSize_fixef4 <- fixef(Sstd2_lognorm_fragSize4)
S_PIE_lognorm_fragSize_fixef4 <- fixef(S_PIE_lognorm_fragSize4)
Scov_lognorm_fragSize_fixef4 <- fixef(Scov_lognorm_fragSize4)
Sn_lognorm_fragSize_fixef4 <- fixef(Sn_lognorm_fragSize4)
Schao_lognorm_fragSize_fixef4 <- fixef(S_chao_lognorm_fragSize4)
Nstd_lognorm_fragSize_fixef4 <- fixef(Nstd_lognorm_fragSize4)

Sstd_posterior <- study_sample_posterior  %>% 
  select(-data) %>% 
  unnest(S_std) %>% 
  mutate(response = 'Sstd',
         Sstd_global_slope = Sstd2_lognorm_fragSize_fixef['c.lfs','Estimate'],
         Sstd_upper_slope = Sstd2_lognorm_fragSize_fixef['c.lfs','Q97.5'],
         Sstd_lower_slope = Sstd2_lognorm_fragSize_fixef['c.lfs','Q2.5']) %>% 
  separate(cbmtt, into = c('continent', 'biome', 'Matrix.category', 'time.since.fragmentation', 'taxa'), 
           sep = '_', remove = F) %>% 
  left_join(meta %>% 
              select(-Case.id, -dataset.name, -reference,
                     -country, -coordinates, -sphere.fragment,
                     -sphere.matrix), 
            by = c('continent', 'biome', 'Matrix.category', 'time.since.fragmentation', 'taxa', 'dataset_label'))
Sstd_posterior$time.since.fragmentation <- factor(Sstd_posterior$time.since.fragmentation,
                                                  levels = c('Recent (less than 20 years)',
                                                             'Intermediate (20-100 years)',
                                                             'long (100+ years)'),
                                                  labels = c('< 20 years',
                                                             '20-100 years)',
                                                             '> 100 years)'))

sstd_study_posterior_time_taxa <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(xmin = Sstd_lower_slope, xmax = Sstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = Sstd_posterior,
                      aes(x = S_std + unique(Sstd_global_slope), 
                          y = time.since.fragmentation,
                          fill = taxa
                      ),
                      scale = 1, alpha = 0.3,
                      linetype = 0) +
  scale_fill_brewer(name = 'Taxa', type = 'qual', palette = 'Accent') +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = '',#Time since fragmentation
       x = ''#expression(paste('Study-level slope'))#,
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(1, 1), 
        legend.justification = c(1, 1),
        legend.background = element_blank())

sstd_study_posterior_time_continent <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(xmin = Sstd_lower_slope, xmax = Sstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = Sstd_posterior,
                      aes(x = S_std + unique(Sstd_global_slope), 
                          y = time.since.fragmentation,
                          fill = continent
                      ),
                      scale = 1, alpha = 0.5,
                      linetype = 0) +
  scale_fill_brewer(name = 'Continent', type = 'qual', palette = 'Dark2') +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = '',
       x = ''#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(1, 1), 
        legend.justification = c(1, 1),
        legend.background = element_blank())

sstd_study_posterior_time_biome <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(xmin = Sstd_lower_slope, xmax = Sstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = Sstd_posterior,
                      aes(x = S_std + unique(Sstd_global_slope), 
                          y = time.since.fragmentation,
                          fill = biome
                      ),
                      scale = 1, alpha = 0.5,
                      linetype = 0) +
  scale_fill_brewer(name = 'Biome', type = 'qual', palette = 'Set2') +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = '',
       x = ''#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(1, 1), 
        legend.justification = c(1, 1),
        legend.background = element_blank())

sstd_study_posterior_time_matrix <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(xmin = Sstd_lower_slope, xmax = Sstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = Sstd_posterior,
                      aes(x = S_std + unique(Sstd_global_slope), 
                          y = time.since.fragmentation,
                          fill = Matrix.category
                      ),
                      scale = 1, alpha = 0.5,
                      linetype = 0) +
  scale_fill_brewer(name = 'Matrix', type = 'qual', palette = 'Set1') +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = '',
       x = ''#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(1, 1), 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

sstd_study_posterior_time <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(xmin = Sstd_lower_slope, xmax = Sstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = Sstd_posterior,
                      aes(x = S_std + unique(Sstd_global_slope), 
                          y = time.since.fragmentation,
                          # fill = taxa
                      ),
                      scale = 1, alpha = 0.5,
                      linetype = 0) +
  # scale_fill_viridis_d(name = 'Taxa') +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = '',
       x = '',#expression(paste('Study-level slope')),
       subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(1, 1), 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

top <- cowplot::plot_grid(NULL, sstd_study_posterior_time, NULL,
                          ncol = 3,
                          rel_widths = c(0.2,0.6,0.2))

bottom <- cowplot::plot_grid(sstd_study_posterior_time_matrix,
                             sstd_study_posterior_time_taxa,
                             sstd_study_posterior_time_biome,
                             sstd_study_posterior_time_continent)

cowplot::plot_grid(top, bottom, nrow = 2,
                   rel_heights = c(0.6,1),
                   rel_widths = c(0.5, 1)) + 
  cowplot::draw_label('Study-level slope', y = 0.009) +
  cowplot::draw_label('Time since fragmentation', x = 0.007, angle = 90)

ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/study_level_slope_time_group.png',
       width = 250,
       height = 250,
       units = 'mm')

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




Sstd_posterior4 <- cbmtt_posterior  %>% 
  select(-data) %>% 
  unnest(S_std) %>% 
  mutate(response = 'Sstd',
         Sstd_global_slope = Sstd2_lognorm_fragSize_fixef4['c.lfs','Estimate'],
         Sstd_upper_slope = Sstd2_lognorm_fragSize_fixef4['c.lfs','Q97.5'],
         Sstd_lower_slope = Sstd2_lognorm_fragSize_fixef4['c.lfs','Q2.5']) %>% 
  separate(cbmtt, into = c('continent', 'biome', 'Matrix.category', 'time.since.fragmentation', 'taxa'), 
           sep = '_', remove = F) %>% 
  left_join(meta %>% 
              select(-dataset_label, -Case.id, -dataset.name, -reference,
                     -country, -coordinates, -sphere.fragment,
                     -sphere.matrix), 
            by = c('continent', 'biome', 'Matrix.category', 'time.since.fragmentation', 'taxa'))

Sstd_posterior4$time.since.fragmentation <- factor(Sstd_posterior4$time.since.fragmentation,
                                                   levels = c('Recent (less than 20 years)',
                                                              'Intermediate (20-100 years)',
                                                              'long (100+ years)'),
                                                   labels = c('< 20 years',
                                                              '20-100 years)',
                                                              '> 100 years)'))

ggplot() +
  facet_grid(. ~ Matrix.category, scale = 'free') +
  geom_rect(data = Sstd_posterior4 %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(xmin = Sstd_lower_slope, xmax = Sstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = Sstd_posterior4 %>% filter(biome=='forest'),
                      aes(x = S_std + unique(Sstd_global_slope), 
                          y = time.since.fragmentation,
                          fill = taxa
                      ),
                      scale = 1, alpha = 0.5,
                      linetype = 0) +
  scale_fill_brewer(name = 'Taxa', type = 'qual', palette = 'Accent') +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = '',
       x = '',#expression(paste('Study-level slope')),
       subtitle = expression(paste('Posterior samples of cbmtt-level ', S[std], ' fragment area slopes (forests only)'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', #c(1, 1), 
        legend.direction = 'horizontal',
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/cbmtt_level_slope_time.png',
       width = 250,
       height = 200,
       units = 'mm')

Spie_posterior4 <- cbmtt_posterior  %>% 
  select(-data) %>% 
  unnest(S_PIE) %>% 
  mutate(response = 'S_PIE',
         Spie_global_slope = S_PIE_lognorm_fragSize_fixef4['c.lfs','Estimate'],
         Spie_upper_slope = S_PIE_lognorm_fragSize_fixef4['c.lfs','Q97.5'],
         Spie_lower_slope = S_PIE_lognorm_fragSize_fixef4['c.lfs','Q2.5']) %>% 
  separate(cbmtt, into = c('continent', 'biome', 'Matrix.category', 'time.since.fragmentation', 'taxa'), 
           sep = '_', remove = F) %>% 
  left_join(meta %>% 
              select(-dataset_label, -Case.id, -dataset.name, -reference,
                     -country, -coordinates, -sphere.fragment,
                     -sphere.matrix), 
            by = c('continent', 'biome', 'Matrix.category', 'time.since.fragmentation', 'taxa'))

Spie_posterior4$time.since.fragmentation <- factor(Spie_posterior4$time.since.fragmentation,
                                                   levels = c('Recent (less than 20 years)',
                                                              'Intermediate (20-100 years)',
                                                              'long (100+ years)'),
                                                   labels = c('< 20 years',
                                                              '20-100 years)',
                                                              '> 100 years)'))

ggplot() +
  facet_grid(. ~ Matrix.category, scale = 'free') +
  geom_rect(data = Spie_posterior4 %>% distinct(Spie_lower_slope, Spie_upper_slope),
            aes(xmin = Spie_lower_slope, xmax = Spie_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = Spie_posterior4 %>% filter(biome=='forest'),
                      aes(x = S_PIE + unique(Spie_global_slope), 
                          y = time.since.fragmentation,
                          fill = taxa
                      ),
                      scale = 1, alpha = 0.5,
                      linetype = 0) +
  scale_fill_brewer(name = 'Taxa', type = 'qual', palette = 'Accent') +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = '',
       x = '',#expression(paste('Study-level slope')),
       subtitle = expression(paste('Posterior samples of cbmtt-level ', S[PIE], ' fragment area slopes (forests only)'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', #c(1, 1), 
        legend.direction = 'horizontal',
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/spie_cbmtt_level_slope_time.png',
       width = 250,
       height = 200,
       units = 'mm')

NStd_posterior4 <- cbmtt_posterior  %>% 
  select(-data) %>% 
  unnest(Nstd) %>% 
  mutate(response = 'Nstd',
         Nstd_global_slope = Nstd_lognorm_fragSize_fixef4['c.lfs','Estimate'],
         Nstd_upper_slope = Nstd_lognorm_fragSize_fixef4['c.lfs','Q97.5'],
         Nstd_lower_slope = Nstd_lognorm_fragSize_fixef4['c.lfs','Q2.5']) %>% 
  separate(cbmtt, into = c('continent', 'biome', 'Matrix.category', 'time.since.fragmentation', 'taxa'), 
           sep = '_', remove = F) %>% 
  left_join(meta %>% 
              select(-dataset_label, -Case.id, -dataset.name, -reference,
                     -country, -coordinates, -sphere.fragment,
                     -sphere.matrix), 
            by = c('continent', 'biome', 'Matrix.category', 'time.since.fragmentation', 'taxa'))

NStd_posterior4$time.since.fragmentation <- factor(NStd_posterior4$time.since.fragmentation,
                                                   levels = c('Recent (less than 20 years)',
                                                              'Intermediate (20-100 years)',
                                                              'long (100+ years)'),
                                                   labels = c('< 20 years',
                                                              '20-100 years)',
                                                              '> 100 years)'))

ggplot() +
  facet_grid(. ~ Matrix.category, scale = 'free') +
  geom_rect(data = NStd_posterior4 %>% distinct(Nstd_lower_slope, Nstd_upper_slope),
            aes(xmin = Nstd_lower_slope, xmax = Nstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = NStd_posterior4 %>% filter(biome=='forest'),
                      aes(x = Nstd + unique(Nstd_global_slope), 
                          y = time.since.fragmentation,
                          fill = taxa
                      ),
                      scale = 1, alpha = 0.5,
                      linetype = 0) +
  scale_fill_brewer(name = 'Taxa', type = 'qual', palette = 'Accent') +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = '',
       x = '',#expression(paste('Study-level slope')),
       subtitle = expression(paste('Posterior samples of cbmtt-level ', N[std], ' fragment area slopes (forests only)'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', #c(1, 1), 
        legend.direction = 'horizontal',
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/Nstd_cbmtt_level_slope_time.png',
       width = 250,
       height = 200,
       units = 'mm')
