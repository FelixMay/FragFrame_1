# plot study-level slopes: most interested in relationship between Std and Jtu for different times since fragmentation

# load packages and paths: 0_init_dirs_load_packages.R
# get the beta-regression coefficients
source(paste0(path2wd, 'r_scripts/5d_beta_frag_coef_wrangle.R'))

# get the fragment area coefficients
source(paste0(path2wd, 'r_scripts/5a_fragSize_coef_wrangle.R'))

# combine study-level estimates from the different analyses (beta vs alpha)
study_slope_coefs <- Jtu_z1i_group_coefs %>% 
  mutate(jtu_intercept = Intercept, 
         jtu_slope = Slope,
         jtu_upper = Slope_upper,
         jtu_lower = Slope_lower,
         jtu_zoi_slope = zoi_Slope,
         jtu_zoi_upper = zoi_Slope_upper,
         jtu_zoi_lower = zoi_Slope_lower) %>% 
  select(dataset_label, jtu_intercept, jtu_slope, jtu_upper, jtu_lower, jtu_zoi_slope, jtu_zoi_upper, jtu_zoi_lower) %>% 
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
    Sstd_lognorm_fragSize_group_coefs %>% 
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


# contribution of turnover and nestedness to dissimilarity: most studies have turnover > nestedness 
study_slope_coefs %>% 
  filter(jtu_intercept > jne_intercept)
study_slope_coefs %>% # 
  filter(rtu_intercept > rne_intercept)

# change in component with fragment size difference: mostly delta_turnover > delta_nestedness,
study_slope_coefs %>% 
  filter(jtu_slope > jne_slope)
study_slope_coefs %>% 
  filter(rtu_slope > rne_slope)

s_jtu_corr <- cor.test(study_slope_coefs$jtu_slope, study_slope_coefs$Sstd)
s_jne_corr <- cor.test(study_slope_coefs$Sstd, study_slope_coefs$jne_slope)

timeLegend <-
  ggplot() +
  # facet_grid(.~climate) +
  geom_point(data = study_slope_coefs,
             aes(y = jtu_slope, x = Sstd, 
                 colour = time.since.fragmentation
             ), 
             # alpha = 0.8
             ) +
  scale_color_manual(name = 'Time since fragmentation',
                     values = c('20-100 years' = '#969696',
                                '< 20 years' = '#cccccc',
                                '> 100 years' = '#636363')) +
  theme_bw() +
  theme(legend.position = 'right',
        legend.direction = 'vertical',
        # legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.text = element_text(size = 5, face = 'plain', hjust = 0),
        legend.title = element_text(size = 6, face = 'plain'),
        legend.margin = margin(),
        legend.box.spacing = unit(c(0,0,0,0), units = 'mm'),
        legend.key.size = unit(2, units = 'mm'),
        plot.margin = unit(c(0,0,0,0), units = 'mm')
  ) +
  guides(colour = guide_legend(title = 'Time since\nfragmentation', label.hjust = 0))

source(paste0(path2wd, 'r_scripts/99_gg_legend.R'))
time_colour_legend <- gg_legend(timeLegend)

beta_turnover_sstd_slope <-
  ggplot() +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(data = study_slope_coefs,
             aes(x = jtu_slope, y = Sstd, 
                 colour = time.since.fragmentation
             ),
             # alpha = 0.8
             ) +
  stat_smooth(data = study_slope_coefs,
              aes(x = jtu_slope, y = Sstd,
                  # colour = time.since.fragmentation
              ),
              method = 'lm', se = F, colour = 'black', size = 0.75
  ) +
  # annotate('text', x = -Inf, y = -0.02, hjust = -0.5, vjust = 0,
  #          label = paste("paste(italic(rho) == " ,
  #                        round(s_jtu_corr$estimate, 2), " (95*'%'~CI: ",
  #                        round(s_jtu_corr$conf.int[1], 2),
  #                        "~\u2013",
  #                        round(s_jtu_corr$conf.int[2], 2),"))"),
  #          parse = T, size = 1.5) +
  scale_color_manual(name = 'Time since fragmentation',
                     values = c('20-100 years' = '#969696',
                                '< 20 years' = '#cccccc',
                                '> 100 years' = '#636363')) +
  labs(x = 'Standardised turnover ~ fragment size slope',
       y = ''
       # y = expression(paste('Study-level ', S[std], 'slope'))
  ) +
  theme_bw() +
  coord_fixed(ratio = 0.95/0.2) +
  theme(legend.position = 'none',# #
        legend.justification = c(0,0),
        legend.box.spacing = unit(0, units = 'mm'),
        legend.background = element_blank(),
        # panel.grid.minor = element_blank(),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        plot.margin = unit(c(0,4,0,0), units = 'mm'))


beta_nestedness_sstd_study <-
  ggplot() +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(data = study_slope_coefs,
             aes(x = jne_slope, y = Sstd, 
                 colour = time.since.fragmentation),
             # alpha = 0.8
             ) +
  stat_smooth(data = study_slope_coefs,
              aes(x = jne_slope, y = Sstd,
                  # colour = time.since.fragmentation
              ),
              method = 'lm', se = F, colour = 'black', size = 0.75
  ) +
  # annotate('text', x = Inf, y = -0.02, hjust = 1.2, vjust = 0,
  #          label = paste("paste(italic(rho) == " ,
  #                        round(s_jne_corr$estimate, 2), " (95*'%'~CI: ",
  #                        round(s_jne_corr$conf.int[1], 2),
  #                        " - ",
  #                        round(s_jne_corr$conf.int[2], 2),"))"),
  #          parse = T, size = 1.5) +
  scale_color_manual(name = 'Time since fragmentation',
                     values = c('20-100 years' = '#969696',
                                '< 20 years' = '#cccccc',
                                '> 100 years' = '#636363')) +
  labs(x = 'Standardised nestedness ~ fragment size slope',
       y = ''
       # y = expression(paste('Study-level biodiversity measure slope'))
  ) +
  theme_bw() +
  coord_fixed(ratio = 1.5/0.2) +
  theme(legend.position = 'none',
        legend.justification = c(1,1),
        # panel.grid.minor = element_blank(),
        # legend.direction = 'horizontal',
        legend.background = element_blank(), 
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        plot.margin = unit(c(0,4,0,0), units = 'mm'))

right = time_colour_legend
left = cowplot::plot_grid(beta_turnover_sstd_slope, 
                          beta_nestedness_sstd_study,
                          align = 'hv',
                          nrow = 1,
                          labels = 'auto',label_size = 8, label_fontface = 'bold') +
  cowplot::draw_label(expression(paste('Standardised richness ~ fragement size slope')),
                      angle = 90,
                      x = 0.015, y = 0.5, size = 6)
cowplot::plot_grid(left, right, nrow = 1, rel_widths = c(1, 0.14))

# two column size for print version
# setwd for saving locally
ggsave(paste0(path2wd, 'main_results/Fig4.pdf'),
       width = 120, height = 55, units = 'mm')
