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
  mutate(jtu_intercept = Intercept, 
         jtu_slope = Slope,
         jtu_upper = Slope_upper,
         jtu_lower = Slope_lower) %>% 
  select(dataset_label, jtu_intercept, jtu_slope, jtu_upper, jtu_lower) %>% 
  left_join(Rtu_z1i_group_coefs %>% 
              mutate(rtu_intercept = Intercept, 
                     rtu_slope = Slope, 
                     rtu_upper = Slope_upper,
                     rtu_lower = Slope_lower) %>% 
              select(dataset_label, rtu_intercept, rtu_slope, rtu_upper, rtu_lower),
            by = 'dataset_label'
            ) %>% 
  left_join(
    Jne_zi_group_coefs %>% 
      mutate(jne_intercept = Intercept, 
             jne_slope = Slope, 
             jne_upper = Slope_upper,
             jne_lower = Slope_lower) %>% 
      select(dataset_label, jne_intercept, jne_slope, jne_upper, jne_lower),
    by = 'dataset_label'
  ) %>% 
  left_join(
    Rne_zi_group_coefs %>% 
      mutate(rne_intercept = Intercept, 
             rne_slope = Slope, 
             rne_upper = Slope_upper,
             rne_lower = Slope_lower) %>% 
      select(dataset_label, rne_intercept, rne_slope, rne_upper, rne_lower),
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
  left_join(
    Nstd_fragSize_group_coefs %>% 
      mutate(Nstd = Slope,
             Nstd_upper = Slope_upper,
             Nstd_lower = Slope_lower) %>% 
      select(dataset_label, Nstd, Nstd_upper, Nstd_lower),
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


# contribution of turnover and nestedness to dissimilarity: most studies have turnover > nestedness (119/123)
study_slope_coefs %>% 
  filter(jtu_intercept > jne_intercept)
study_slope_coefs %>% 
  filter(rtu_intercept > rne_intercept)

# change in component with fragment size difference: 42/123 have delta_turnover > delta_nestedness, 
# increases to 95/123 for abundance based distance
study_slope_coefs %>% 
  filter(jtu_slope > jne_slope)
study_slope_coefs %>% 
  filter(rtu_slope > rne_slope)

s_jtu_corr <- cor.test(study_slope_coefs$Sstd, study_slope_coefs$jtu_slope)
s_jne_corr <- cor.test(study_slope_coefs$Sstd, study_slope_coefs$jne_slope)

timeLegend <-
  ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(y = jtu_slope, x = Sstd, 
                 colour = time.since.fragmentation
             )) +
  scale_color_manual(name = 'Time since fragmentation', 
                     values = c('20-100 years' = '#6996b3',
                                '< 20 years' = '#c1e7ff',
                                '> 100 years' = '#004c6d')) +
  theme_bw() +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        # legend.justification = c(1,1),
        legend.background = element_blank())

source('~/Dropbox/1current/R_random/functions/gg_legend.R')
time_colour_legend <- gg_legend(timeLegend)

beta_turnover_sstd_slope <-
ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(x = jtu_slope, y = Sstd, 
                  colour = time.since.fragmentation
                 )) +
  stat_smooth(data = study_slope_coefs,
            aes(x = jtu_slope, y = Sstd
                ),
            method = 'lm', se = F, colour = 'black'
            ) +
  # annotate('text', x = Inf, y = -0.02, hjust = 1.025, vjust = 0,
  #          label = paste("paste(italic(rho) == " , 
  #                        round(s_jtu_corr$estimate, 2), " (95*'%'~CI: ",
  #                        round(s_jtu_corr$conf.int[1], 2),
  #                        " ~`???` ",
  #                        round(s_jtu_corr$conf.int[2], 2),"))"),
  #          parse = T, size = 2.5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(name = 'Time since fragmentation', 
                     values = c('20-100 years' = '#6996b3',
                                '< 20 years' = '#c1e7ff',
                                '> 100 years' = '#004c6d')) +
  labs(x = 'Study-level turnover slope',
       y = ''
       # y = expression(paste('Study-level ', S[std], 'slope'))
       ) +
  theme_bw() +
  theme(legend.position = 'none')


beta_nestedness_sstd_study <-
ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
           aes(x = jne_slope, y = Sstd, 
               colour = time.since.fragmentation)) +
  stat_smooth(data = study_slope_coefs,
              aes(x = jne_slope, y = Sstd,
                  # colour = time.since.fragmentation
                  ),
              method = 'lm', se = F, colour = 'black'
              ) +
  # annotate('text', x = Inf, y = -0.02, hjust = 1.2, vjust = 0,
  #          label = paste("paste(italic(rho) == " , 
  #                        round(s_jne_corr$estimate, 2), " (95*'%'~CI: ",
  #                        round(s_jne_corr$conf.int[1], 2),
  #                        " ~`???`~",
  #                        round(s_jne_corr$conf.int[2], 2),"))"),
  #          parse = T, size = 2.5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(values = c('20-100 years' = '#6996b3',
                                '< 20 years' = '#c1e7ff',
                                '> 100 years' = '#004c6d')) +
  labs(x = 'Study-level nestedness slope',
       y = ''
       # y = expression(paste('Study-level biodiversity measure slope'))
       ) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.justification = c(1,1),
        legend.direction = 'horizontal',
        legend.background = element_blank())

top = time_colour_legend
bottom = cowplot::plot_grid(beta_turnover_sstd_slope, 
                   beta_nestedness_sstd_study,
                   align = 'hv',
                   nrow = 1,
                   labels = 'auto') +
  cowplot::draw_label(expression(paste('Study-level ', S[std], 'slope')),
                      angle = 90,
                      x = 0.02)

cowplot::plot_grid(top, bottom, nrow = 2, rel_heights = c(0.1, 1))
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/Fig4_colour.png', width = 170, height = 80, units = 'mm')
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/fig2_taxa_colour.png', width = 250, height = 80, units = 'mm')

# relationship between the intercepts of turnover and nestedness tell us about their
# relative contribution to total dissimilarity
cor.test(study_slope_coefs$jtu_intercept, study_slope_coefs$jne_intercept)
with(study_slope_coefs, summary(lm(jne_intercept ~ jtu_intercept)))
# beta_turnover_nestedness_intercept <-
  ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(x = jtu_intercept, y = jne_intercept, 
                 # colour = time.since.fragmentation
             )) +
  stat_smooth(data = study_slope_coefs,
              aes(x = jtu_intercept, y = jne_intercept
              ),
              method = 'lm', se = F, colour = 'black'
  ) +
  # annotate('text', x = Inf, y = -0.02, hjust = 1.025, vjust = 0,
  #          label = paste("paste(italic(rho) == " , 
  #                        round(s_jtu_corr$estimate, 2), " (95*'%'~CI: ",
  #                        round(s_jtu_corr$conf.int[1], 2),
  #                        " ~`???` ",
  #                        round(s_jtu_corr$conf.int[2], 2),"))"),
  #          parse = T, size = 2.5) +
  # geom_hline(yintercept = 0, lty = 2) +
  # geom_vline(xintercept = 0, lty = 2) +
  # scale_color_manual(name = 'Time since fragmentation', 
  #                    values = c('20-100 years' = '#6996b3',
  #                               '< 20 years' = '#c1e7ff',
  #                               '> 100 years' = '#004c6d')) +
  # labs(x = 'Study-level turnover slope',
  #      y = ''
  #      # y = expression(paste('Study-level ', S[std], 'slope'))
  # ) +
  theme_bw() +
  theme(legend.position = 'none')

# beta_turnover_sstd_slope_colour <-
  ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(x = jtu_slope, y = Sstd, 
                 colour = 'Sstd'
             )) +
  # all alpha metrics?
  geom_point(data = study_slope_coefs,
             aes(x = jtu_slope, y = Nstd, 
                 colour = 'Nstd'
             )) +
  geom_point(data = study_slope_coefs,
             aes(x = jtu_slope, y = Spie, 
                 colour = 'S_PIE'
             )) +
    geom_point(data = study_slope_coefs,
               aes(x = jtu_slope, y = Schao, 
                   colour = 'Schao'
               )) +
    geom_point(data = study_slope_coefs,
               aes(x = jtu_slope, y = Sn, 
                   colour = 'Sn'
               )) +
    geom_point(data = study_slope_coefs,
               aes(x = jtu_slope, y = Scov, 
                   colour = 'Scov'
               )) +  
  stat_smooth(data = study_slope_coefs,
              aes(x = jtu_slope, y = Sstd, 
                  colour = 'Sstd'
              ),
              method = 'lm', se = F
  ) +
    stat_smooth(data = study_slope_coefs,
                aes(x = jtu_slope, y = Nstd, 
                    colour = 'Nstd'
                ),
                method = 'lm', se = F
    ) +
    stat_smooth(data = study_slope_coefs,
                aes(x = jtu_slope, y = Spie, 
                    colour = 'S_PIE'
                ),
                method = 'lm', se = F
    ) +  
  # annotate('text', x = Inf, y = -0.02, hjust = 1, vjust = 0,
  #          label = paste("paste(italic(rho) == " , 
  #                        round(s_jtu_corr$estimate, 2), " (95*'%'~CI: ",
  #                        round(s_jtu_corr$conf.int[1], 2),
  #                        " ~`???` ",
  #                        round(s_jtu_corr$conf.int[2], 2),"))"),
  #          parse = T, size = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  scale_colour_manual(name= '', 
                        values = c('Nstd' = '#1b9e77', 
                                   'Sstd' = '#d95f02', 
                                   'S_PIE' = '#7570b3', 
                                   'Scov' = '#a6761d', 
                                   'Schao' = '#e6ab02',
                                   'Sn' = '#e7298a')) +
  labs(x = 'Study-level turnover slope',
       y = ''
       # y = expression(paste('Study-level ', S[std], 'slope'))
  ) +
  theme_bw() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank())


beta_nestedness_sstd_study <-
  ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(x = jne_slope, y = Sstd)) +
  stat_smooth(data = study_slope_coefs,
              aes(x = jne_slope, y = Sstd),
              method = 'lm', se = F) +
  annotate('text', x = Inf, y = -0.02, hjust = 1.2, vjust = 0,
           label = paste("paste(italic(rho) == " , 
                         round(s_jne_corr$estimate, 2), " (95*'%'~CI: ",
                         round(s_jne_corr$conf.int[1], 2),
                         " ~`???`~",
                         round(s_jne_corr$conf.int[2], 2),"))"),
           parse = T, size = 2) +
  
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = 'Study-level nestedness slope',
       y = ''
       # y = expression(paste('Study-level biodiversity measure slope'))
  ) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.justification = c(1,1),
        legend.direction = 'horizontal',
        legend.background = element_blank())

cowplot::plot_grid(beta_turnover_sstd_slope, 
                   beta_nestedness_sstd_study,
                   align = 'hv',
                   nrow = 1,
                   labels = 'auto') +
  cowplot::draw_label(expression(paste('Study-level ', S[std], 'slope')),
                      angle = 90, x = 0.02)

ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/Fig5.png', width = 170, height = 80, units = 'mm')
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/fig2_taxa_colour.png', width = 250, height = 80, units = 'mm')


study_slope_coefs %>% 
  group_by(time.since.fragmentation) %>% 
  summarise(sum(continent8=='Europe')
            )
