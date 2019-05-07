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
                                                             '20-100 years)',
                                                             '> 100 years)'))

Sstd_posterior$Matrix.category <- factor(Sstd_posterior$Matrix.category,
                                                  levels = c('light filter', 'intermediate', 'harsh filter'),
                                                  labels = c('Light',
                                                             'Intermediate',
                                                             'Harsh'))
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
  geom_density_ridges_gradient(data = Sstd_posterior,
                      aes(x = S_std + unique(Sstd_global_slope),
                          y = time.since.fragmentation,
                          fill = stat(quantile)
                      ),
                      quantiles = c(0.25, 0.75),
                      calc_ecdf = T,
                      scale = 0.95, alpha = 0.5,
                      linetype = 0) +
  # stat_density_ridges(data = Sstd_posterior,
  #                     geom = 'density_ridges_gradient',
  #                     aes(x = S_std + unique(Sstd_global_slope), 
  #                     y = time.since.fragmentation,
  #                     fill = 0.5 - abs(0.5 -..ecdf..)),
  #                     calc_ecdf = T,
  #                     scale = 0.95, alpha = 0.5) +
  geom_point(data = Sstd_posterior,
             aes(x = S_std + unique(Sstd_global_slope), 
                 y = time.since.fragmentation),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 3) +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(xmin = Sstd_lower_slope, xmax = Sstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = '',
       x = '',#expression(paste('Study-level slope')),
       subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc',
                               '#969696',
                               '#cccccc'),
                    labels = c('< 25%', 
                               '25 - 75%',
                               '> 75%')) +
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


sstd_study_posterior_matrx <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = Sstd_posterior,
                               aes(x = S_std + unique(Sstd_global_slope),
                                   y = Matrix.category,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.25, 0.75),
                               calc_ecdf = T,
                               scale = 0.95, alpha = 0.5,
                               linetype = 0) +
  # stat_density_ridges(data = Sstd_posterior,
  #                     geom = 'density_ridges_gradient',
  #                     aes(x = S_std + unique(Sstd_global_slope), 
  #                     y = time.since.fragmentation,
  #                     fill = 0.5 - abs(0.5 -..ecdf..)),
  #                     calc_ecdf = T,
  #                     scale = 0.95, alpha = 0.5) +
  geom_point(data = Sstd_posterior,
             aes(x = S_std + unique(Sstd_global_slope), 
                 y = Matrix.category),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 3) +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(xmin = Sstd_lower_slope, xmax = Sstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Matrix filter',
       x = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc',
                               '#969696',
                               '#cccccc'),
                    labels = c('< 25%', 
                               '25 - 75%',
                               '> 75%')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(1, 1), 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+


sstd_study_posterior_biome <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = Sstd_posterior,
                               aes(x = S_std + unique(Sstd_global_slope),
                                   y = biome,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.25, 0.75),
                               calc_ecdf = T,
                               scale = 0.95, alpha = 0.5,
                               linetype = 0) +
  # stat_density_ridges(data = Sstd_posterior,
  #                     geom = 'density_ridges_gradient',
  #                     aes(x = S_std + unique(Sstd_global_slope), 
  #                     y = time.since.fragmentation,
  #                     fill = 0.5 - abs(0.5 -..ecdf..)),
  #                     calc_ecdf = T,
  #                     scale = 0.95, alpha = 0.5) +
  geom_point(data = Sstd_posterior,
             aes(x = S_std + unique(Sstd_global_slope), 
                 y = biome),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 3) +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(xmin = Sstd_lower_slope, xmax = Sstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Biome',
       x = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc',
                               '#969696',
                               '#cccccc'),
                    labels = c('< 25%', 
                               '25 - 75%',
                               '> 75%')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(1, 1), 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

sstd_study_posterior_taxa <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = Sstd_posterior,
                               aes(x = S_std + unique(Sstd_global_slope),
                                   y = taxa,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.25, 0.75),
                               calc_ecdf = T,
                               scale = 0.95, alpha = 0.5,
                               linetype = 0) +
  # stat_density_ridges(data = Sstd_posterior,
  #                     geom = 'density_ridges_gradient',
  #                     aes(x = S_std + unique(Sstd_global_slope), 
  #                     y = time.since.fragmentation,
  #                     fill = 0.5 - abs(0.5 -..ecdf..)),
  #                     calc_ecdf = T,
  #                     scale = 0.95, alpha = 0.5) +
  geom_point(data = Sstd_posterior,
             aes(x = S_std + unique(Sstd_global_slope), 
                 y = taxa),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 3) +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(xmin = Sstd_lower_slope, xmax = Sstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Taxa',
       x = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc',
                               '#969696',
                               '#cccccc'),
                    labels = c('< 25%', 
                               '25 - 75%',
                               '> 75%')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(1, 1), 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

sstd_study_posterior_continent <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = Sstd_posterior,
                               aes(x = S_std + unique(Sstd_global_slope),
                                   y = continent8,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.25, 0.75),
                               calc_ecdf = T,
                               scale = 0.95, alpha = 0.5,
                               linetype = 0) +
  # stat_density_ridges(data = Sstd_posterior,
  #                     geom = 'density_ridges_gradient',
  #                     aes(x = S_std + unique(Sstd_global_slope), 
  #                     y = time.since.fragmentation,
  #                     fill = 0.5 - abs(0.5 -..ecdf..)),
  #                     calc_ecdf = T,
  #                     scale = 0.95, alpha = 0.5) +
  geom_point(data = Sstd_posterior,
             aes(x = S_std + unique(Sstd_global_slope), 
                 y = continent8),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 3) +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(xmin = Sstd_lower_slope, xmax = Sstd_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_vline(data = Sstd_posterior,
             aes(xintercept = Sstd_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Continent',
       x = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc',
                               '#969696',
                               '#cccccc'),
                    labels = c('< 25%', 
                               '25 - 75%',
                               '> 75%')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(1, 1), 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

cowplot::plot_grid(sstd_study_posterior_time, 
                   sstd_study_posterior_matrx, 
                   sstd_study_posterior_taxa,
                   sstd_study_posterior_biome,
                   sstd_study_posterior_continent)
