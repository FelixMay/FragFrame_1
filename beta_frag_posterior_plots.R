# get study-level posterior samples
library(tidyverse)
library(brms)
library(ggridges)

load('~/Dropbox/1current/fragmentation_synthesis/results/Jtu_z1i_frag_beta-5098650.Rdata')
load('~/Dropbox/1current/fragmentation_synthesis/results/Rtu_z1i_frag_beta-5098692.Rdata')

load('~/Dropbox/1current/fragmentation_synthesis/results/Jne_zi_fragSize.Rdata')
load('~/Dropbox/1current/fragmentation_synthesis/results/Rne_zi_fragSize.Rdata')

frag_beta <- read_csv('~/Dropbox/Frag Database (new)/files_datapaper/Analysis/2_betapart_frag_fcont_10_mabund_as_is.csv')

# get the metadata
meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag_beta <- frag_beta %>% 
  group_by(dataset_label, sample_design, method, frag_x) %>% 
  mutate(pair_group = paste0(frag_x, '_g')) %>% 
  ungroup() %>% 
  # centre covariate before fitting
  mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))

# frag_beta <- left_join(frag_beta, 
#                        meta, by = 'dataset_label')

# study-levels (use model with fewest missing values)
study_levels <- Jtu_z1i_fragSize$data %>% 
  as_tibble() %>% 
  distinct(dataset_label) %>% 
  mutate(level = dataset_label) %>%
  nest(level)

study_sample_posterior <- study_levels %>%
  mutate(Jtu_z1i = purrr::map(data, ~posterior_samples(Jtu_z1i_fragSize, 
                                                     pars = paste('r_dataset_label[', as.character(.x$level), ',cl10ra]', sep=''),
                                                     exact = TRUE,
                                                     subset = floor(runif(n = 1000,
                                                                          min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         Jne_zi = purrr::map(data, ~posterior_samples(Jne_zi_fragSize, 
                                                    pars = paste('r_dataset_label[', as.character(.x$level), ',cl10ra]', sep=''),
                                                    exact = TRUE,
                                                    subset = floor(runif(n = 1000,
                                                                         min = 1, max = 2000))) %>% unlist() %>% as.numeric()))


Jtu_z1i_fixef <- fixef(Jtu_z1i_fragSize)
Jne_zi_fixef <- fixef(Jne_zi_fragSize)

Jtu_z1i_post <- study_sample_posterior  %>% 
  select(-data) %>% 
  unnest(Jtu_z1i) %>% 
  mutate(response = 'Turnover',
         Jtu_global_slope = Jtu_z1i_fixef['cl10ra','Estimate'],
         Jtu_upper_slope = Jtu_z1i_fixef['cl10ra','Q97.5'],
         Jtu_lower_slope = Jtu_z1i_fixef['cl10ra','Q2.5']) %>% 
  left_join(meta, 
            by = 'dataset_label')

Jne_zi_post <- study_sample_posterior  %>% 
  select(-data) %>% 
  unnest(Jne_zi) %>% 
  mutate(response = 'Nestedness',
         Jne_global_slope = Jne_zi_fixef['cl10ra','Estimate'],
         Jne_upper_slope = Jne_zi_fixef['cl10ra','Q97.5'],
         Jne_lower_slope = Jne_zi_fixef['cl10ra','Q2.5']) %>% 
  left_join(meta, 
            by = 'dataset_label')

Jtu_z1i_post$time.since.fragmentation <- factor(Jtu_z1i_post$time.since.fragmentation,
                                                levels = c('Recent (less than 20 years)',
                                                           'Intermediate (20-100 years)',
                                                           'long (100+ years)'),
                                                labels = c('Recent (less than 20 years)',
                                                           'Intermediate (20-100 years)',
                                                           'Long (100+ years)'))

Jne_zi_post$time.since.fragmentation <- factor(Jne_zi_post$time.since.fragmentation,
                                                  levels = c('Recent (less than 20 years)',
                                                             'Intermediate (20-100 years)',
                                                             'long (100+ years)'),
                                                  labels = c('Recent (less than 20 years)',
                                                             'Intermediate (20-100 years)',
                                                             'Long (100+ years)'))


ggplot() +
  facet_grid( ~ biome, scale = 'free') +
  geom_rect(data = Jtu_z1i_post %>% distinct(Jtu_lower_slope, Jtu_upper_slope),
            aes(xmin = Jtu_lower_slope, xmax = Jtu_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = Jtu_z1i_post,
                      aes(x = Jtu_z1i + unique(Jtu_global_slope), 
                          y = interaction(Matrix.category, time.since.fragmentation),
                          fill = taxa
                      ),
                      scale = 1, alpha = 0.6,
                      linetype = 0) +
  # scale_fill_viridis_d(name = 'Taxa') +
  geom_vline(data = Jtu_z1i_post,
             aes(xintercept = Jtu_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Time since fragmentation & matrix category',
       x = expression(paste('Study-level ', J[tu], ' change'))#,
       # tag = 'A'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.5, 0.2)) #+

ggplot() +
  facet_grid(biome ~ taxa, scales = 'free') +
  geom_rect(data = Jtu_z1i_post %>% distinct(Jtu_lower_slope, Jtu_upper_slope),
            aes(xmin = Jtu_lower_slope, xmax = Jtu_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = Jtu_z1i_post,
                      aes(x = Jtu_z1i + unique(Jtu_global_slope), 
                          y = time.since.fragmentation,
                          fill = Matrix.category
                      ),
                      scale = 0.9, alpha = 0.6,
                      linetype = 0) +
  # scale_linetype_manual(name = 'Realm', values = c('Marine' = 0, 'Terrestrial' = 1, 'Freshwater' = 3)) +
  # scale_fill_manual(name = 'Taxa', values = taxa_col) +
  geom_vline(data = Jtu_z1i_post,
             aes(xintercept = Jtu_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Time since fragmentation',
       x = expression(paste('Study-level turnover slope'))#,
       # tag = 'A'
  ) +
  # scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.1, 0.1),
        legend.direction = 'horizontal',
        #legend.background = element_blank()
        ) #+

ggplot() +
  facet_grid(biome ~ taxa, scales = 'free') +
  geom_rect(data = Jne_zi_post %>% distinct(Jne_lower_slope, Jne_upper_slope),
            aes(xmin = Jne_lower_slope, xmax = Jne_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.3) +
  geom_density_ridges(data = Jne_zi_post,
                      aes(x = Jne_zi + unique(Jne_global_slope), 
                          y = time.since.fragmentation,
                          fill = Matrix.category
                      ),
                      scale = 0.9, alpha = 0.6,
                      linetype = 0) +
  # scale_linetype_manual(name = 'Realm', values = c('Marine' = 0, 'Terrestrial' = 1, 'Freshwater' = 3)) +
  # scale_fill_manual(name = 'Taxa', values = taxa_col) +
  geom_vline(data = Jne_zi_post,
             aes(xintercept = Jne_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Time since fragmentation',
       x = expression(paste('Study-level nestedness slope'))#,
       # tag = 'A'
  ) +
  # scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.1, 0.1),
        legend.direction = 'horizontal',
        #legend.background = element_blank()
  ) #+
