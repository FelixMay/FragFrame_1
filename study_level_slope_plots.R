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

beta_turnover_sstd_slope <-
ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(x = jtu_slope, y = Sstd, 
                 # colour = time.since.fragmentation
                 )) +
  stat_smooth(data = study_slope_coefs,
            aes(x = jtu_slope, y = Sstd#, colour = time.since.fragmentation
                ),
            method = 'lm', se = F
            ) +
  annotate('text', x = -Inf, y = Inf, hjust = -0.2, vjust = 1.4,
           label = paste("paste(italic(rho) == " , 
                         round(s_jtu_corr$estimate, 2), " (95*'%'~CI: ",
                         round(s_jtu_corr$conf.int[1], 2),
                         " ~`—` ",
                         round(s_jtu_corr$conf.int[2], 2),"))"),
           parse = T, size = 5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
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
  annotate('text', x = -Inf, y = Inf, hjust = -0.2, vjust = 1.4,
           label = paste("paste(italic(rho) == " , 
                         round(s_jne_corr$estimate, 2), " (95*'%'~CI: ",
                         round(s_jne_corr$conf.int[1], 2),
                         " ~`—` ",
                         round(s_jne_corr$conf.int[2], 2),"))"),
           parse = T, size = 5) +
  
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
                   nrow = 2) +
  cowplot::draw_label(expression(paste('Study-level ', S[std], 'slope')),
                      angle = 90, x = 0.01)

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/Fig5.png',
       # width = 110, height = 180, units = 'mm')


# beta_turnover2 <- 
ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(x = rtu_slope, y = jtu_slope),
             size = 0.5) +
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

# beta_nestedness2 <- 
ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(x = rne_slope, y = jne_slope),
             size = 0.5) +
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

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/SStd_beta_components_study_level.png', 
#        width = 250, height = 150, units = 'mm')


ggplot() +
  # facet_grid(taxa~continent8) +
  geom_point(data = study_slope_coefs,
             aes(x = jtu_slope, y = jne_slope),#x = jtu_slope, y = jne_slope
             size = 0.5) +
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

# study-level alpha components: 
S_N <- lm(Sstd ~ Nstd, data = study_slope_coefs)
S_N_corr <- cor.test(study_slope_coefs$Nstd[!is.na(study_slope_coefs$Sstd)], 
                        study_slope_coefs$Sstd[!is.na(study_slope_coefs$Sstd)])
S_Spie <- lm(Sstd ~ Spie, data = study_slope_coefs)
S_Spie_corr <- cor.test(study_slope_coefs$Spie[!is.na(study_slope_coefs$Sstd)], 
                   study_slope_coefs$Sstd[!is.na(study_slope_coefs$Sstd)])

# add indicator for difference from zero (within studies)
study_slope_coefs <- study_slope_coefs %>% 
  mutate(effect = ifelse((Sstd_lower > 0 & Sstd_upper > 0) & (Spie_lower > 0 & Spie_upper > 0) & (Nstd_lower > 0 & Nstd_upper > 0),
                             'S, S_PIE & N decay',
                             ifelse((Sstd_lower > 0 & Sstd_upper > 0) & (Spie_lower < 0 & Spie_upper > 0) & (Nstd_lower < 0 & Nstd_upper > 0),
                                    'S decay',
                                    ifelse((Sstd_lower < 0 & Sstd_upper > 0) & (Spie_lower < 0 & Spie_upper > 0) & (Nstd_lower > 0 & Nstd_upper > 0),
                                           'N decay', 
                                           ifelse((Sstd_lower < 0 & Sstd_upper > 0) & (Spie_lower > 0 & Spie_upper > 0) & (Nstd_lower < 0 & Nstd_upper > 0),
                                                  'S_PIE decay',
                                                  ifelse((Sstd_lower > 0 & Sstd_upper > 0) & (Spie_lower > 0 & Spie_upper > 0) & (Nstd_lower < 0 & Nstd_upper > 0),
                                                         'S & S_PIE decay',
                                                         ifelse((Sstd_lower > 0 & Sstd_upper > 0) & (Spie_lower < 0 & Spie_upper > 0) & (Nstd_lower > 0 & Nstd_upper > 0),
                                                                'S & N decay', 'random')))))))

# order for plotting
study_slope_coefs$effect <- factor(study_slope_coefs$effect,
                                       levels = c('random', 'S decay', 'S_PIE decay', 'N decay', 
                                                  'S & N decay', 'S & S_PIE decay', 'S, S_PIE & N decay'))

alpha_study1 <- ggplot() +
  # facet_grid(taxa~continent8) +
  geom_point(data = study_slope_coefs %>% filter(!is.na(Sstd)),
             aes(x = Nstd, y = Sstd, colour = effect),
             size = 1.5) +
  stat_smooth(data = study_slope_coefs,
              aes(x = Nstd, y = Sstd),
              method = 'lm', se = F) +
  annotate('text', x = -Inf, y = Inf, hjust = -1, vjust = 1.4,
           label = paste("paste(italic(rho) == " , 
                         round(S_N_corr$estimate, 2), ")"),  
           parse = T, size = 5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  labs(x = expression(paste('Study-level ',N[std],' slope')),
       y = expression(paste('Study-level ',S[std],' slope'))) +
  scale_colour_manual(name= '', 
                      values = c('random' = '#1b9e77', 
                                 'S decay' = '#d95f02', 
                                 'N decay' = '#7570b3', 
                                 'S_PIE decay' = '#a6761d', 
                                 'S & S_PIE decay' = '#e6ab02',
                                 'S & N decay' = '#e7298a',
                                 'S, S_PIE & N decay' = '#666666')) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.justification = c(0,0),
        legend.direction = 'horizontal',
        legend.background = element_blank()) +
  guides(colour = guide_legend(nrow = 7))


alpha_study2 <- ggplot() +
  # facet_grid(taxa~continent8) +
  geom_point(data = study_slope_coefs %>% filter(!is.na(Sstd)),
             aes(x = Spie, y = Sstd, colour = effect),
             size = 1.5) +
  stat_smooth(data = study_slope_coefs,
              aes(x = Spie, y = Sstd),
              method = 'lm', se = F) +
  annotate('text', x = -Inf, y = Inf, hjust = -1, vjust = 1.4,
           label = paste("paste(italic(rho) == " , 
                         round(S_Spie_corr$estimate, 2), ")"),  
           parse = T, size = 5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  labs(x = expression(paste('Study-level ',S[PIE],' slope')),
       y = '' #expression(paste('Study-level ',S[std],' slope'))
       ) +
  scale_colour_manual(name= '', 
                      values = c('random' = '#1b9e77', 
                                 'S decay' = '#d95f02', 
                                 'N decay' = '#7570b3', 
                                 'S_PIE decay' = '#a6761d', 
                                 'S & S_PIE decay' = '#e6ab02',
                                 'S & N decay' = '#e7298a',
                                 'S, S_PIE & N decay' = '#666666')) +
  theme_bw() +
  theme(legend.position = c(0, 0.5),
        legend.justification = c(0,0),
        legend.direction = 'horizontal',
        legend.background = element_blank()) +
  guides(colour = guide_legend(nrow = 7))

cowplot::plot_grid(alpha_study1, alpha_study2,
                   nrow = 1, align = 'hv')

ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/FigSx_alpha_study_level_reln.png',
       width = 290, height = 150, units = 'mm')

# turnover_alpha_slopes <- 
  ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(x = jtu_slope, y = Sstd),
             size = 0.5, colour = 'blue') +
  geom_point(data = study_slope_coefs,
             aes(x = jtu_slope, y = Spie),
             size = 0.5, colour = 'black') +
  geom_point(data = study_slope_coefs,
               aes(x = jtu_slope, y = Nstd),
               size = 0.5, colour = 'dark red') +
  stat_smooth(data = study_slope_coefs,
            aes(x = jtu_slope, y = Sstd),
            method = 'gam', se = F) +
  stat_smooth(data = study_slope_coefs,
              aes(x = jtu_slope, y = Spie),
              method = 'gam', se = F, colour = 'black') +
    stat_smooth(data = study_slope_coefs,
                aes(x = jtu_slope, y = Nstd),
                method = 'gam', se = F, colour = 'red') +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = 'Study-level turnover slope',
       y = expression(paste('Study-level ', S[std], ', ', S[PIE],
                            ', ', N[std], ' slope'))) +
  theme_bw() +
  theme(legend.position = c(0, 0),
        legend.justification = c(0,0),
        legend.direction = 'horizontal',
        legend.background = element_blank())

# nestedness_alpha_slopes <- 
  ggplot() +
    # facet_grid(.~climate) +
    geom_point(data = study_slope_coefs,
               aes(x = jne_slope, y = Sstd),
               size = 0.5, colour = 'blue') +
    geom_point(data = study_slope_coefs,
               aes(x = jne_slope, y = Spie),
               size = 0.5, colour = 'black') +
    geom_point(data = study_slope_coefs,
               aes(x = jne_slope, y = Nstd),
               size = 0.5, colour = 'dark red') +
    stat_smooth(data = study_slope_coefs,
                aes(x = jne_slope, y = Sstd),
                method = 'gam', se = F) +
    stat_smooth(data = study_slope_coefs,
                aes(x = jne_slope, y = Spie),
                method = 'gam', se = F, colour = 'black') +
    stat_smooth(data = study_slope_coefs,
                aes(x = jne_slope, y = Nstd),
                method = 'gam', se = F, colour = 'red') +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = 'Study-level nestedness slope',
         y = expression(paste('Study-level ', S[std], ', ', S[PIE],
                              ', ', N[std], ' slope'))) +
    theme_bw() +
    theme(legend.position = c(0, 0),
          legend.justification = c(0,0),
          legend.direction = 'horizontal',
          legend.background = element_blank())
  
  ggplot() +
    # facet_grid(.~climate) +
    geom_point(data = study_slope_coefs,
               aes(x = rtu_slope, y = Sstd),
               size = 0.5, colour = 'blue') +
    geom_point(data = study_slope_coefs,
               aes(x = rtu_slope, y = Spie),
               size = 0.5, colour = 'black') +
    geom_point(data = study_slope_coefs,
               aes(x = rtu_slope, y = Nstd),
               size = 0.5, colour = 'dark red') +
    stat_smooth(data = study_slope_coefs,
                aes(x = rtu_slope, y = Sstd),
                method = 'gam', se = F) +
    stat_smooth(data = study_slope_coefs,
                aes(x = rtu_slope, y = Spie),
                method = 'gam', se = F, colour = 'black') +
    stat_smooth(data = study_slope_coefs,
                aes(x = rtu_slope, y = Nstd),
                method = 'gam', se = F, colour = 'red') +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = 'Study-level balanced abundance slope',
         y = expression(paste('Study-level ', S[std], ', ', S[PIE],
                              ', ', N[std], ' slope'))) +
    theme_bw() +
    theme(legend.position = c(0, 0),
          legend.justification = c(0,0),
          legend.direction = 'horizontal',
          legend.background = element_blank())
  
  # turnover_alpha_slopes <- 
  ggplot() +
    # facet_grid(.~climate) +
    geom_point(data = study_slope_coefs,
               aes(x = rne_slope, y = Sstd),
               size = 0.5, colour = 'blue') +
    geom_point(data = study_slope_coefs,
               aes(x = rne_slope, y = Spie),
               size = 0.5, colour = 'black') +
    geom_point(data = study_slope_coefs,
               aes(x = rne_slope, y = Nstd),
               size = 0.5, colour = 'dark red') +
    stat_smooth(data = study_slope_coefs,
                aes(x = rne_slope, y = Sstd),
                method = 'gam', se = F) +
    stat_smooth(data = study_slope_coefs,
                aes(x = rne_slope, y = Spie),
                method = 'gam', se = F, colour = 'black') +
    stat_smooth(data = study_slope_coefs,
                aes(x = rne_slope, y = Nstd),
                method = 'gam', se = F, colour = 'dark red') +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = 'Study-level abundance gradient slope',
         y = expression(paste('Study-level ', S[std], ', ', S[PIE],
                              ', ', N[std], ' slope'))) +
    theme_bw() +
    theme(legend.position = c(0, 0),
          legend.justification = c(0,0),
          legend.direction = 'horizontal',
          legend.background = element_blank())


  
  # alternate df for different plots
long_df <- bind_rows(
    Jtu_z1i_group_coefs %>% 
    select(Intercept, Intercept_lower, Intercept_upper, dataset_label, Slope, Slope_lower, Slope_upper) %>% 
    mutate(response = 'jtu'),
    Rtu_z1i_group_coefs %>% 
      select(Intercept, Intercept_lower, Intercept_upper, dataset_label, Slope, Slope_lower, Slope_upper) %>% 
      mutate(response = 'rtu'),
    Jne_zi_group_coefs %>% 
      select(Intercept, Intercept_lower, Intercept_upper, dataset_label, Slope, Slope_lower, Slope_upper) %>% 
      mutate(response = 'jne'),
    Rne_zi_group_coefs %>% 
      select(Intercept, Intercept_lower, Intercept_upper, dataset_label, Slope, Slope_lower, Slope_upper) %>% 
      mutate(response = 'rne'),
    Sstd2_lognorm_fragSize_group_coefs %>% 
      select(Intercept, Intercept_lower, Intercept_upper, dataset_label, Slope, Slope_lower, Slope_upper) %>% 
      mutate(response = 'Sstd'),
    Nstd_fragSize_group_coefs %>% 
      select(Intercept, Intercept_lower, Intercept_upper, dataset_label, Slope, Slope_lower, Slope_upper) %>% 
      mutate(response = 'Nstd'),
    S_PIE_fragSize_group_coefs %>% 
      select(Intercept, Intercept_lower, Intercept_upper, dataset_label, Slope, Slope_lower, Slope_upper) %>% 
      mutate(response = 'Spie'),
    Sn_lognorm_fragSize_group_coefs %>% 
      select(Intercept, Intercept_lower, Intercept_upper, dataset_label, Slope, Slope_lower, Slope_upper) %>% 
      mutate(response = 'Sn'),
    Schao_lognorm_fragSize_group_coefs %>% 
    select(Intercept, Intercept_lower, Intercept_upper, dataset_label, Slope, Slope_lower, Slope_upper) %>% 
    mutate(response = 'Schao'),
    Scov_lognorm_fragSize_group_coefs %>% 
    select(Intercept, Intercept_lower, Intercept_upper, dataset_label, Slope, Slope_lower, Slope_upper) %>% 
    mutate(response = 'Scov')) %>% 
    mutate(effect = ifelse((Slope_lower < 0 & Slope_upper > 0), 'random',
                           ifelse((Slope_lower > 0 & Slope_upper > 0), 'decay', 'down'))) %>% 
    left_join(meta, 
              by = 'dataset_label') %>% 
  group_by(response) %>% 
  mutate(n_random = sum(effect=='random'),
         n_decay = sum(effect=='decay')) %>% 
  ungroup()

dissimilarity = c('jtu', 'rtu', 'jne', 'rne')

long_df$time.since.fragmentation <- factor(long_df$time.since.fragmentation,
                                                     levels = c('Recent (less than 20 years)',
                                                                'Intermediate (20-100 years)',
                                                                'long (100+ years)'),
                                                     labels = c('< 20 years',
                                                                '20-100 years',
                                                                '> 100 years'))

long_df$Matrix.category <- factor(long_df$Matrix.category,
                                            levels = c('light filter', 'intermediate', 'harsh filter'),
                                            labels = c('Light', 'Intermediate', 'Harsh'))

long_df$biome <- factor(long_df$biome,
                                  levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                                  labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))
long_df$taxa <- factor(long_df$taxa,
                                 levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                                 labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

global_long <- bind_rows(
  Jtu_z1i_fixef %>% 
  as_tibble() %>% 
  mutate(term = rownames(Jtu_z1i_fixef),
         response = 'jtu'),
  Rtu_z1i_fixef %>% 
    as_tibble() %>% 
    mutate(term = rownames(Rtu_z1i_fixef),
           response = 'rtu'),
  Jne_zi_fixef %>% 
    as_tibble() %>% 
    mutate(term = rownames(Jne_zi_fixef),
           response = 'jne'),
  Rne_zi_fixef %>% 
    as_tibble() %>% 
    mutate(term = rownames(Rne_zi_fixef),
           response = 'rne'),
  Nstd_lognorm_fragSize_fixef %>% 
    as_tibble() %>% 
    mutate(term = rownames(Nstd_lognorm_fragSize_fixef),
           response = 'Nstd'),
  Sn_lognorm_fragSize_fixef %>% 
    as_tibble() %>% 
    mutate(term = rownames(Sn_lognorm_fragSize_fixef),
           response = 'Sn'),
  Sstd2_lognorm_fragSize_fixef %>% 
    as_tibble() %>% 
    mutate(term = rownames(Sstd2_lognorm_fragSize_fixef),
           response = 'Sstd'),
  S_PIE_lognorm_fragSize_fixef %>% 
    as_tibble() %>% 
    mutate(term = rownames(S_PIE_lognorm_fragSize_fixef),
           response = 'Spie'),
  Schao_lognorm_fragSize_fixef %>% 
    as_tibble() %>% 
    mutate(term = rownames(Schao_lognorm_fragSize_fixef),
           response = 'Schao'),
  Scov_lognorm_fragSize_fixef %>% 
    as_tibble() %>% 
    mutate(term = rownames(Scov_lognorm_fragSize_fixef),
           response = 'Scov'))
  

long_df %>% 
    filter(response%in%dissimilarity) %>% 
    ggplot() +
    facet_wrap(~response) +
    geom_density(aes(x = Slope)) +
    geom_vline(xintercept = 0, lty = 2)
  
long_df %>% 
  filter(!response%in%dissimilarity) %>% 
  ggplot() +
  facet_wrap(~response) +
  geom_histogram(aes(x = Slope, fill = effect)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_vline(data = global_long %>% filter(!response%in%dissimilarity & term=='c.lfs'), 
             aes(xintercept = Estimate)) +
  geom_rect(data = global_long %>% filter(!response%in%dissimilarity & term=='c.lfs'), 
            aes(ymin = -Inf, ymax = Inf, xmin = Q2.5, xmax = Q97.5),
            alpha = 0.3) +
  # add proportion of random and decay
  geom_text(data = long_df %>% filter(!response%in%dissimilarity) %>% 
             distinct(response, n_random, n_decay, effect),
             aes(x = Inf, y = Inf, hjust = 1.8, vjust = 1.4,
                 group = response, colour = 'random',
                 label = paste("n[random] == ", #[Frag.~size]
                         n_random)),  
           parse = T) +
  geom_text(data = long_df %>% filter(!response%in%dissimilarity) %>% 
              distinct(response, n_random, n_decay, effect),
            aes(x = Inf, y = Inf, hjust = 2, vjust = 2.4,
                group = response, colour = 'decay',
                label = paste("n[decay] == ", #[Frag.~size]
                              n_decay)),  
            parse = T) +
  scale_fill_manual(values = c('random' = '#1b9e77', 'decay' = '#7570b3')) +
  scale_colour_manual(values = c('random' = '#1b9e77', 'decay' = '#7570b3')) +
  labs(x = 'Slope estimate',
       y = 'Number of studies') +
  theme_bw() +
  theme(legend.position = 'none')

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/FigSx_study_level_slope_distr.png',
#        width = 290, height = 200, units = 'mm')

long_df %>% 
  filter(!response%in%dissimilarity & !is.na(time.since.fragmentation)) %>% 
  ggplot() +
  facet_wrap(~response) +
  geom_histogram(aes(x = biome, fill = effect), stat = 'count')

long_df %>% 
  filter(!response%in%dissimilarity & !is.na(time.since.fragmentation)) %>% 
  ggplot() +
  facet_grid(effect~response) +
  geom_histogram(aes(x = time.since.fragmentation, fill = biome), stat = 'count')

# want the maximum difference in area for each dataset_label
long_df %>% 
  filter(!response%in%dissimilarity) %>% 
  left_join(frag %>% 
  group_by(dataset_label) %>% 
  summarise(max_diff = max(frag_size_num) - min(frag_size_num),
            ratio = max(frag_size_num)/ min(frag_size_num)) %>% 
    ungroup(),
  by= 'dataset_label') %>% 
  ggplot() +
  facet_grid(~response) +
  geom_point(aes(x = ratio, y = Slope, colour = effect)) +
  scale_x_continuous(trans = 'log10')


y = Jtu_z1i_fS$data$repl
x = Jtu_z1i_fS$data$cl10ra
group = Jtu_z1i_fS$data$dataset_label
yrep = posterior_predict(Jtu_z1i_fS, summary = F)
bayesplot::ppc_error_scatter_avg_vs_x(y, yrep, x)
bayesplot::ppc_error_hist(y, yrep[1:10,])
