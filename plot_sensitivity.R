# code to plot coefficient estimates of models fit to data standardised in different ways
library(tidyverse)

# code to wrangle the coefs ready to inspect
source('~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/fragSize_coef_wrangle_4_sensitivity.R')

fixed_effects <- bind_cols(
  # reference estimates
  Nstd_ref_fixef %>% 
    as_tibble() %>% 
    mutate(Nstd_ref = Estimate[2],
           Nstd_ref_lower = Q2.5[2],
           Nstd_ref_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Nstd_ref, Nstd_ref_lower, Nstd_ref_upper),
  Sstd2_ref_fixef %>% 
    as_tibble() %>% 
    mutate(Sstd_ref = Estimate[2],
           Sstd_ref_lower = Q2.5[2],
           Sstd_ref_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Sstd_ref, Sstd_ref_lower, Sstd_ref_upper),
  S_PIE_ref_fixef %>% 
    as_tibble() %>% 
    mutate(S_PIE_ref = Estimate[2],
           S_PIE_ref_lower = Q2.5[2],
           S_PIE_ref_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(S_PIE_ref, S_PIE_ref_lower, S_PIE_ref_upper),
  # case 1
  Nstd_sens1_fixef %>% 
    as_tibble() %>% 
    mutate(Nstd_sens1 = Estimate[2],
           Nstd_sens1_lower = Q2.5[2],
           Nstd_sens1_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Nstd_sens1, Nstd_sens1_lower, Nstd_sens1_upper),
  Sstd2_sens1_fixef %>% 
    as_tibble() %>% 
    mutate(Sstd_sens1 = Estimate[2],
           Sstd_sens1_lower = Q2.5[2],
           Sstd_sens1_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Sstd_sens1, Sstd_sens1_lower, Sstd_sens1_upper),
  S_PIE_sens1_fixef %>% 
    as_tibble() %>% 
    mutate(S_PIE_sens1 = Estimate[2],
           S_PIE_sens1_lower = Q2.5[2],
           S_PIE_sens1_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(S_PIE_sens1, S_PIE_sens1_lower, S_PIE_sens1_upper),
  # case 3  
  Nstd_sens3_fixef %>% 
    as_tibble() %>% 
    mutate(Nstd_sens3 = Estimate[2],
           Nstd_sens3_lower = Q2.5[2],
           Nstd_sens3_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Nstd_sens3, Nstd_sens3_lower, Nstd_sens3_upper),
  Sstd2_sens3_fixef %>% 
    as_tibble() %>% 
    mutate(Sstd_sens3 = Estimate[2],
           Sstd_sens3_lower = Q2.5[2],
           Sstd_sens3_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Sstd_sens3, Sstd_sens3_lower, Sstd_sens3_upper),
  S_PIE_sens3_fixef %>% 
    as_tibble() %>% 
    mutate(S_PIE_sens3 = Estimate[2],
           S_PIE_sens3_lower = Q2.5[2],
           S_PIE_sens3_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(S_PIE_sens3, S_PIE_sens3_lower, S_PIE_sens3_upper),
  # case 8
  Nstd_sens8_fixef %>% 
    as_tibble() %>% 
    mutate(Nstd_sens8 = Estimate[2],
           Nstd_sens8_lower = Q2.5[2],
           Nstd_sens8_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Nstd_sens8, Nstd_sens8_lower, Nstd_sens8_upper),
  Sstd2_sens8_fixef %>% 
    as_tibble() %>% 
    mutate(Sstd_sens8 = Estimate[2],
           Sstd_sens8_lower = Q2.5[2],
           Sstd_sens8_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Sstd_sens8, Sstd_sens8_lower, Sstd_sens8_upper),
  S_PIE_sens8_fixef %>% 
    as_tibble() %>% 
    mutate(S_PIE_sens8 = Estimate[2],
           S_PIE_sens8_lower = Q2.5[2],
           S_PIE_sens8_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(S_PIE_sens8, S_PIE_sens8_lower, S_PIE_sens8_upper),
  # case 11
  Nstd_sens11_fixef %>% 
    as_tibble() %>% 
    mutate(Nstd_sens11 = Estimate[2],
           Nstd_sens11_lower = Q2.5[2],
           Nstd_sens11_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Nstd_sens11, Nstd_sens11_lower, Nstd_sens11_upper),
  Sstd2_sens11_fixef %>% 
    as_tibble() %>% 
    mutate(Sstd_sens11 = Estimate[2],
           Sstd_sens11_lower = Q2.5[2],
           Sstd_sens11_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Sstd_sens11, Sstd_sens11_lower, Sstd_sens11_upper),
  S_PIE_sens11_fixef %>% 
    as_tibble() %>% 
    mutate(S_PIE_sens11 = Estimate[2],
           S_PIE_sens11_lower = Q2.5[2],
           S_PIE_sens11_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(S_PIE_sens11, S_PIE_sens11_lower, S_PIE_sens11_upper))
                           

# study level
study_coefs <- left_join(
  # reference fits first
  Nstd_ref_grp_coefs %>% 
    mutate(Nstd_ref = Slope,
           Nstd_ref_upper = Slope_upper, 
           Nstd_ref_lower = Slope_upper) %>% 
    select(dataset_label, Nstd_ref, Nstd_ref_lower, Nstd_ref_upper),
  Sstd2_ref_grp_coefs %>% 
    mutate(Sstd_ref = Slope,
           Sstd_ref_upper = Slope_upper, 
           Sstd_ref_lower = Slope_upper) %>% 
    select(dataset_label, Sstd_ref, Sstd_ref_lower, Sstd_ref_upper),
  by = 'dataset_label') %>% 
  inner_join(
  S_PIE_ref_grp_coefs %>% 
    mutate(S_PIE_ref = Slope,
           S_PIE_ref_upper = Slope_upper, 
           S_PIE_ref_lower = Slope_upper) %>% 
    select(dataset_label, S_PIE_ref, S_PIE_ref_lower, S_PIE_ref_upper),
  by = 'dataset_label') %>% 
  inner_join(
    # case 1
    Nstd_sens1_grp_coefs %>% 
      mutate(Nstd_sens1 = Slope,
             Nstd_sens1_upper = Slope_upper, 
             Nstd_sens1_lower = Slope_upper) %>% 
      select(dataset_label, Nstd_sens1, Nstd_sens1_lower, Nstd_sens1_upper),
    by = 'dataset_label') %>% 
    left_join(
      Sstd2_sens1_grp_coefs %>% 
        mutate(Sstd_sens1 = Slope,
               Sstd_sens1_upper = Slope_upper, 
               Sstd_sens1_lower = Slope_upper) %>% 
        select(dataset_label, Sstd_sens1, Sstd_sens1_lower, Sstd_sens1_upper),  
      by = 'dataset_label') %>% 
    inner_join(
      S_PIE_sens1_grp_coefs %>% 
        mutate(S_PIE_sens1 = Slope,
               S_PIE_sens1_upper = Slope_upper, 
               S_PIE_sens1_lower = Slope_upper) %>% 
        select(dataset_label, S_PIE_sens1, S_PIE_sens1_lower, S_PIE_sens1_upper),
      by = 'dataset_label') %>% 
    inner_join(
      # case 3
      Nstd_sens3_grp_coefs %>% 
        mutate(Nstd_sens3 = Slope,
               Nstd_sens3_upper = Slope_upper, 
               Nstd_sens3_lower = Slope_upper) %>% 
        select(dataset_label, Nstd_sens3, Nstd_sens3_lower, Nstd_sens3_upper),
      by = 'dataset_label') %>% 
    left_join(
      Sstd2_sens3_grp_coefs %>% 
        mutate(Sstd_sens3 = Slope,
               Sstd_sens3_upper = Slope_upper, 
               Sstd_sens3_lower = Slope_upper) %>% 
        select(dataset_label, Sstd_sens3, Sstd_sens3_lower, Sstd_sens3_upper),  
      by = 'dataset_label') %>% 
    inner_join(
      S_PIE_sens3_grp_coefs %>% 
        mutate(S_PIE_sens3 = Slope,
               S_PIE_sens3_upper = Slope_upper, 
               S_PIE_sens3_lower = Slope_upper) %>% 
        select(dataset_label, S_PIE_sens3, S_PIE_sens3_lower, S_PIE_sens3_upper),
      by = 'dataset_label') %>% 
    inner_join(
      # case 8
      Nstd_sens8_grp_coefs %>% 
        mutate(Nstd_sens8 = Slope,
               Nstd_sens8_upper = Slope_upper, 
               Nstd_sens8_lower = Slope_upper) %>% 
        select(dataset_label, Nstd_sens8, Nstd_sens8_lower, Nstd_sens8_upper),
      by = 'dataset_label') %>% 
    left_join(
      Sstd2_sens8_grp_coefs %>% 
        mutate(Sstd_sens8 = Slope,
               Sstd_sens8_upper = Slope_upper, 
               Sstd_sens8_lower = Slope_upper) %>% 
        select(dataset_label, Sstd_sens8, Sstd_sens8_lower, Sstd_sens8_upper),  
      by = 'dataset_label') %>% 
    inner_join(
      S_PIE_sens8_grp_coefs %>% 
        mutate(S_PIE_sens8 = Slope,
               S_PIE_sens8_upper = Slope_upper, 
               S_PIE_sens8_lower = Slope_upper) %>% 
        select(dataset_label, S_PIE_sens8, S_PIE_sens8_lower, S_PIE_sens8_upper),
      by = 'dataset_label') %>% 
    inner_join(
      # case 11
      Nstd_sens11_grp_coefs %>% 
        mutate(Nstd_sens11 = Slope,
               Nstd_sens11_upper = Slope_upper, 
               Nstd_sens11_lower = Slope_upper) %>% 
        select(dataset_label, Nstd_sens11, Nstd_sens11_lower, Nstd_sens11_upper),
      by = 'dataset_label') %>% 
    left_join(
      Sstd2_sens11_grp_coefs %>% 
        mutate(Sstd_sens11 = Slope,
               Sstd_sens11_upper = Slope_upper, 
               Sstd_sens11_lower = Slope_upper) %>% 
        select(dataset_label, Sstd_sens11, Sstd_sens11_lower, Sstd_sens11_upper),  
      by = 'dataset_label') %>% 
    inner_join(
      S_PIE_sens11_grp_coefs %>% 
        mutate(S_PIE_sens11 = Slope,
               S_PIE_sens11_upper = Slope_upper, 
               S_PIE_sens11_lower = Slope_upper) %>% 
        select(dataset_label, S_PIE_sens11, S_PIE_sens11_lower, S_PIE_sens11_upper),
      by = 'dataset_label')

# %>% 
#   mutate(rho_estimate = cor.test(jtu_norm, jtu_beta, method = 'spearman')$estimate %>% signif(digits = 2),
#          p_value = cor.test(jtu_norm, jtu_beta, method = 'spearman')$p.value %>% signif(digits = 2))


Sstd_robust <- ggplot() +
  # reference versus case 1
  geom_point(data = study_coefs,
             aes(x = Sstd_ref, y = Sstd_sens1, colour = 'sens1'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = Sstd_ref, y = Sstd_sens1, colour = 'sens1'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = Sstd_ref_lower, xmax = Sstd_ref_upper, y = Sstd_sens1, colour = 'sens1'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = Sstd_ref, ymin = Sstd_sens1_lower, ymax = Sstd_sens1_upper, colour = 'sens1')) +
  # reference versus case 3
  geom_point(data = study_coefs,
             aes(x = Sstd_ref, y = Sstd_sens3, colour = 'sens3'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = Sstd_ref, y = Sstd_sens3, colour = 'sens3'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = Sstd_ref_lower, xmax = Sstd_ref_upper, y = Sstd_sens3, colour = 'sens3'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = Sstd_ref, ymin = Sstd_sens3_lower, ymax = Sstd_sens3_upper, colour = 'sens3')) +
  # reference versus case 8
  geom_point(data = study_coefs,
             aes(x = Sstd_ref, y = Sstd_sens8, colour = 'sens8'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = Sstd_ref, y = Sstd_sens8, colour = 'sens8'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = Sstd_ref_lower, xmax = Sstd_ref_upper, y = Sstd_sens8, colour = 'sens8'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = Sstd_ref, ymin = Sstd_sens8_lower, ymax = Sstd_sens8_upper, colour = 'sens8')
                 ) +
  # reference versus case 11
  geom_point(data = study_coefs,
             aes(x = Sstd_ref, y = Sstd_sens11, colour = 'sens11'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = Sstd_ref, y = Sstd_sens11, colour = 'sens11'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = Sstd_ref_lower, xmax = Sstd_ref_upper, y = Sstd_sens11, colour = 'sens11'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = Sstd_ref, ymin = Sstd_sens11_lower, ymax = Sstd_sens11_upper, colour = 'sens11')) +
  # 1:1 line and zero lines
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  # geom_vline(xintercept = 0, lty = 2) +
  # geom_hline(yintercept = 0, lty = 2) +
  scale_colour_manual(name = 'Standardisation',
                      guide = F,
                      values = c('sens1' = '#003f5c',
                                 'sens3' = '#7a5195',
                                 'sens8' = '#ef5675',
                                 'sens11' = '#ffa600'),
                      labels = c('1', '2', '3', '4')) +
  labs(x = '',
       y = '',
       subtitle = expression(paste('Species richness'))) +# (', S[std], ')
  theme_bw()


Nstd_robust <-
  ggplot() +
  # reference versus case 1
  geom_point(data = study_coefs,
             aes(x = Nstd_ref, y = Nstd_sens1, colour = 'sens1'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = Nstd_ref, y = Nstd_sens1, colour = 'sens11'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = Nstd_ref_lower, xmax = Nstd_ref_upper, y = Nstd_sens1, colour = 'sens1'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = Nstd_ref, ymin = Nstd_sens1_lower, ymax = Nstd_sens1_upper, colour = 'sens1')) +
  # reference versus case 3
  geom_point(data = study_coefs,
             aes(x = Nstd_ref, y = Nstd_sens3, colour = 'sens3'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = Nstd_ref, y = Nstd_sens3, colour = 'sens3'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = Nstd_ref_lower, xmax = Nstd_ref_upper, y = Nstd_sens3, colour = 'sens3'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = Nstd_ref, ymin = Nstd_sens3_lower, ymax = Nstd_sens3_upper, colour = 'sens3')) +
  # reference versus case 8
  geom_point(data = study_coefs,
             aes(x = Nstd_ref, y = Nstd_sens8, colour = 'sens8'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = Nstd_ref, y = Nstd_sens8, colour = 'sens8'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = Nstd_ref_lower, xmax = Nstd_ref_upper, y = Nstd_sens8, colour = 'sens8'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = Nstd_ref, ymin = Nstd_sens8_lower, ymax = Nstd_sens8_upper, colour = 'sens8')) +
  # reference versus case 11
  geom_point(data = study_coefs,
             aes(x = Nstd_ref, y = Nstd_sens11, colour = 'sens11'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = Nstd_ref, y = Nstd_sens11, colour = 'sens11'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = Nstd_ref_lower, xmax = Nstd_ref_upper, y = Nstd_sens11, colour = 'sens11'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = Nstd_ref, ymin = Nstd_sens11_lower, ymax = Nstd_sens11_upper, colour = 'sens11')) +
  # 1:1 line and zero lines
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  # geom_vline(xintercept = 0, lty = 2) +
  # geom_hline(yintercept = 0, lty = 2) +
  scale_colour_manual(name = 'Standardisation',
                      values = c('sens1' = '#003f5c',
                                 'sens3' = '#7a5195',
                                 'sens8' = '#ef5675',
                                 'sens11' = '#ffa600'),
                      labels = c('1', '2', '3', '4')) +
  labs(x = '',
       y = '',
       subtitle = expression(paste('Number of individuals'))) +# (', N[std], ')
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.background = element_blank(),
        # legend.box = element_blank(),
        legend.justification = c(1,0)) +
  guides(colour = guide_legend(nrow = 2))
  

S_PIE_robust <- ggplot() +
  # reference versus case 1
  geom_point(data = study_coefs,
             aes(x = S_PIE_ref, y = S_PIE_sens1, colour = 'sens1'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = S_PIE_ref, y = S_PIE_sens1, colour = 'sens1'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = S_PIE_ref_lower, xmax = S_PIE_ref_upper, y = S_PIE_sens1, colour = 'sens1'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = S_PIE_ref, ymin = S_PIE_sens1_lower, ymax = S_PIE_sens1_upper, colour = 'sens1')) +
  # reference versus case 3
  geom_point(data = study_coefs,
             aes(x = S_PIE_ref, y = S_PIE_sens3, colour = 'sens3'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = S_PIE_ref, y = S_PIE_sens3, colour = 'sens3'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = S_PIE_ref_lower, xmax = S_PIE_ref_upper, y = S_PIE_sens3, colour = 'sens3'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = S_PIE_ref, ymin = S_PIE_sens3_lower, ymax = S_PIE_sens3_upper, colour = 'sens3')) +
  # reference versus case 8
  geom_point(data = study_coefs,
             aes(x = S_PIE_ref, y = S_PIE_sens8, colour = 'sens8'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = S_PIE_ref, y = S_PIE_sens8, colour = 'sens8'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = S_PIE_ref_lower, xmax = S_PIE_ref_upper, y = S_PIE_sens8, colour = 'sens8'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = S_PIE_ref, ymin = S_PIE_sens8_lower, ymax = S_PIE_sens8_upper, colour = 'sens8')) +
  # reference versus case 11
  geom_point(data = study_coefs,
             aes(x = S_PIE_ref, y = S_PIE_sens11, colour = 'sens11'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = S_PIE_ref, y = S_PIE_sens11, colour = 'sens11'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = S_PIE_ref_lower, xmax = S_PIE_ref_upper, y = S_PIE_sens11, colour = 'sens11'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = S_PIE_ref, ymin = S_PIE_sens11_lower, ymax = S_PIE_sens11_upper, colour = 'sens11')) +
  # 1:1 line and zero lines
  scale_colour_manual(name = 'Standardisation',
                      guide = F,
                      values = c('sens1' = '#003f5c',
                                 'sens3' = '#7a5195',
                                 'sens8' = '#ef5675',
                                 'sens11' = '#ffa600'),
                      labels = c('1', '2', '3', '4')) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  # geom_vline(xintercept = 0, lty = 2) +
  # geom_hline(yintercept = 0, lty = 2) +
  labs(x = '',
       y = '',
       subtitle = expression(paste('Evenness'))) +# (', S[PIE], ')
  theme_bw()

cowplot::plot_grid(Nstd_robust,
                   Sstd_robust,
                   S_PIE_robust,
                   nrow = 1, align = 'hv',
                   labels = 'auto') +
  cowplot::draw_label('Slope estimate for reference standardisation',
                      y = 0.03, size = 13) +
  cowplot::draw_label('Slope estimate for\nalternate standardisation',
                      x = 0.02, angle = 90, size = 13)

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/figSx_robust_main.png',
#        width = 250,
#        height = 80,
#        units = 'mm')
