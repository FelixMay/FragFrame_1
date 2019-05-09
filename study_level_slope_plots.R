# plot study-level slopes: most interested in relationship between Std and Jtu for different times since fragmentation

rm(list=ls())
library(tidyverse)
library(brms)
# get the beta-regression coefficients
source('~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/beta_frag_coef_wrangle.R')

# get the fragment area coefficients
source('~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/fragSize_coef_wrangle.R')

# combine study-level estimates from the different analyses (beta vs alpha)
study_slope_coefs <- Jtu_z1i_group_coefs %>% 
  mutate(jtu_slope = Slope,
         jtu_upper = Slope_upper,
         jtu_lower = Slope_lower) %>% 
  select(dataset_label, jtu_slope, jtu_upper, jtu_lower) %>% 
  left_join(Rtu_z1i_group_coefs %>% 
              mutate(rtu_slope = Slope, 
                     rtu_upper = Slope_upper,
                     rtu_lower = Slope_lower) %>% 
              select(dataset_label, rtu_slope, rtu_upper, rtu_lower),
            by = 'dataset_label'
            ) %>% 
  left_join(
    Jne_zi_group_coefs %>% 
      mutate(jne_slope = Slope, 
             jne_upper = Slope_upper,
             jne_lower = Slope_lower) %>% 
      select(dataset_label, jne_slope, jne_upper, jne_lower),
    by = 'dataset_label'
  ) %>% 
  left_join(
    Rne_zi_group_coefs %>% 
      mutate(rne_slope = Slope, 
             rne_upper = Slope_upper,
             rne_lower = Slope_lower) %>% 
      select(dataset_label, rne_slope, rne_upper, rne_lower),
    by = 'dataset_label'
  ) %>% 
  left_join(
    Sstd2_lognorm_fragSize_group_coefs %>% 
      mutate(Sstd = Slope,
             Sstd_upper = Slope_upper,
             Sstd_lower = Slope_lower) %>% 
      select(dataset_label, Sstd, Sstd_upper, Sstd_lower),
    by = 'dataset_label'
  ) %>% 
  left_join(
    S_PIE_fragSize_group_coefs %>% 
      mutate(Spie = Slope,
             Spie_upper = Slope_upper,
             Spie_lower = Slope_lower) %>% 
      select(dataset_label, Spie, Spie_upper, Spie_lower),
    by = 'dataset_label'
  ) %>% 
  left_join(meta, by = 'dataset_label')

study_slope_coefs$time.since.fragmentation <- factor(study_slope_coefs$time.since.fragmentation,
                                                 levels = c('Recent (less than 20 years)',
                                                            'Intermediate (20-100 years)',
                                                            'long (100+ years)'),
                                                 labels = c('< 20 years',
                                                            '20-100 years',
                                                            '> 100 years'))

study_slope_coefs$Matrix.category <- factor(study_slope_coefs$Matrix.category,
                                        levels = c('light filter', 'intermediate', 'harsh filter'),
                                        labels = c('Light', 'Intermediate', 'Harsh'))

study_slope_coefs$biome <- factor(study_slope_coefs$biome,
                              levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                              labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))
study_slope_coefs$taxa <- factor(study_slope_coefs$taxa,
                             levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                             labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

study_slope_coefs <- study_slope_coefs %>% 
  unite(col = 'frag_matrix', c(sphere.fragment, sphere.matrix),
        remove = F, sep = ', ')


beta_turnover_sstd_slope <- ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(x = jtu_slope, y = Sstd),
             size = 2) +
  # geom_linerange(data = study_slope_coefs,
  #                aes(x = jtu_slope, ymin = Sstd_lower, ymax = Sstd_upper,
  #                    colour = time.since.fragmentation),
  #                lwd = 0.1, alpha = 0.5) +
  # geom_errorbarh(data = study_slope_coefs,
  #                aes(xmin = jtu_lower, xmax = jtu_upper, y = Sstd,
  #                    colour = time.since.fragmentation),
  #                lwd = 0.1, height = 0, alpha = 0.5) +
  # geom_point(data = study_slope_coefs,
  #            aes(x = rtu_slope, y = Sstd, colour = time.since.fragmentation),
  #            shape = 2) +
stat_smooth(data = study_slope_coefs,
            aes(x = jtu_slope, y = Sstd),
            method = 'gam', se = F) +
  # stat_smooth(data = study_slope_coefs,
  #             aes(x = jtu_slope, y = Sstd, colour = continent8),
  #             method = 'lm', se = F,
  #             linetype = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = 'Study-level turnover slope',
       y = expression(paste('Study-level ', S[std], ' slope'))) +
  theme_bw() +
  theme(legend.position = c(0, 0),
        legend.justification = c(0,0),
        legend.direction = 'horizontal',
        legend.background = element_blank())

beta_nestedness_sstd_study <- ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(x = jne_slope, y = Sstd),
             size = 2) +
  # geom_linerange(data = study_slope_coefs,
  #                aes(x = jtu_slope, ymin = Sstd_lower, ymax = Sstd_upper,
  #                    colour = time.since.fragmentation),
  #                lwd = 0.1, alpha = 0.5) +
  # geom_errorbarh(data = study_slope_coefs,
  #                aes(xmin = jtu_lower, xmax = jtu_upper, y = Sstd,
  #                    colour = time.since.fragmentation),
  #                lwd = 0.1, height = 0, alpha = 0.5) +
  # geom_point(data = study_slope_coefs,
  #            aes(x = rtu_slope, y = Sstd, colour = time.since.fragmentation),
  #            shape = 2) +
  stat_smooth(data = study_slope_coefs,
              aes(x = jne_slope, y = Sstd),
              method = 'gam', se = F) +
  # stat_smooth(data = study_slope_coefs,
  #             aes(x = jne_slope, y = Sstd, colour = continent8),
  #             method = 'lm', se = F,
  #             linetype = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = 'Study-level nestedness component slope',
       y = expression(paste('Study-level ', S[std], ' slope'))) +
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.direction = 'horizontal',
        legend.background = element_blank())

cowplot::plot_grid(beta_turnover_sstd_slope, 
                   beta_nestedness_sstd_study,
                   nrow = 2)

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/SStd_beta_components_study_level.png', 
#        width = 250, height = 150, units = 'mm')


beta_turnover2 <- ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(x = rtu_slope, y = jtu_slope),
             size = 2) +
  # geom_linerange(data = study_slope_coefs,
  #                aes(x = jtu_slope, ymin = Sstd_lower, ymax = Sstd_upper,
  #                    colour = time.since.fragmentation),
  #                lwd = 0.1, alpha = 0.5) +
  # geom_errorbarh(data = study_slope_coefs,
  #                aes(xmin = jtu_lower, xmax = jtu_upper, y = Sstd,
  #                    colour = time.since.fragmentation),
  #                lwd = 0.1, height = 0, alpha = 0.5) +
  # geom_point(data = study_slope_coefs,
  #            aes(x = rtu_slope, y = Sstd, colour = time.since.fragmentation),
  #            shape = 2) +
stat_smooth(data = study_slope_coefs,
            aes(x = rtu_slope, y = jtu_slope),
            method = 'gam', se = F) +
  # stat_smooth(data = study_slope_coefs,
  #             aes(x = rtu_slope, y = jtu_slope, colour = continent8),
  #             method = 'lm', se = F,
  #             linetype = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = expression(paste('Study-level balanced abundance slope')),
       y = 'Study-level turnover slope') +
  theme_bw() +
  theme(legend.position = c(0, 0),
        legend.justification = c(0,0),
        legend.direction = 'horizontal',
        legend.background = element_blank())

beta_nestedness2 <- ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(x = rne_slope, y = jne_slope),
             size = 2) +
  # geom_linerange(data = study_slope_coefs,
  #                aes(x = jtu_slope, ymin = Sstd_lower, ymax = Sstd_upper,
  #                    colour = time.since.fragmentation),
  #                lwd = 0.1, alpha = 0.5) +
  # geom_errorbarh(data = study_slope_coefs,
  #                aes(xmin = jtu_lower, xmax = jtu_upper, y = Sstd,
  #                    colour = time.since.fragmentation),
  #                lwd = 0.1, height = 0, alpha = 0.5) +
  # geom_point(data = study_slope_coefs,
  #            aes(x = rtu_slope, y = Sstd, colour = time.since.fragmentation),
  #            shape = 2) +
stat_smooth(data = study_slope_coefs,
            aes(x = rne_slope, y = jne_slope),
            method = 'gam', se = F) +
  # stat_smooth(data = study_slope_coefs,
  #             aes(x = rtu_slope, y = jtu_slope, colour = continent8),
  #             method = 'lm', se = F,
  #             linetype = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = expression(paste('Study-level abundance gradient slope')),
       y = 'Study-level nestedness slope') +
  theme_bw() +
  theme(legend.position = c(0, 0),
        legend.justification = c(0,0),
        legend.direction = 'horizontal',
        legend.background = element_blank())

cowplot::plot_grid(beta_turnover_sstd_slope, 
                   beta_nestedness_sstd_study,
beta_turnover2, 
                   beta_nestedness2,
                   nrow = 2)

ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/SStd_beta_components_study_level.png', 
       width = 250, height = 150, units = 'mm')


ggplot() +
  facet_grid(taxa~continent8) +
  geom_point(data = study_slope_coefs,
             aes(x = jtu_slope, y = jne_slope),
             size = 2) +
  stat_smooth(data = study_slope_coefs,
              aes(x = jtu_slope, y = jne_slope),
              method = 'gam', se = F) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = expression(paste('Study-level turnover slope')),
       y = 'Study-level nestedness slope') +
  theme_bw() +
  theme(legend.position = c(0, 0),
        legend.justification = c(0,0),
        legend.direction = 'horizontal',
        legend.background = element_blank())

ggplot() +
  # facet_grid(taxa~continent8) +
  geom_point(data = study_slope_coefs,
             aes(x = Spie, y = Sstd),
             size = 2) +
  stat_smooth(data = study_slope_coefs,
              aes(x = Spie, y = Sstd),
              method = 'gam', se = F) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = expression(paste('Study-level ',S[PIE],' slope')),
       y = expression(paste('Study-level ',S[std],' slope'))) +
  theme_bw() +
  theme(legend.position = c(0, 0),
        legend.justification = c(0,0),
        legend.direction = 'horizontal',
        legend.background = element_blank())
