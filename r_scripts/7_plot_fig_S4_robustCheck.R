# code to plot coefficient estimates of models fit to data where the size of large fragments 
# was imputed differently (i.e., different multiplicative factors: 2 or 100, 10 used for main text)
# Or where non-integer data were handled differently (rounded) as opposed to left as is in main text

# load packages and paths: 0_init_dirs_load_packages.R

# code to wrangle the coefs ready to inspect
## NB: this code does not run without user having model fits locally (i.e.,
# they are not in the intermediate_results folder)
source(paste0(path2wd, 'r_scripts/5b_fragSize_coef_wrangle_4_robust_results.R'))

fixed_effects <- bind_cols(
  # reference estimates
  Nstd_ref_fixef %>% 
    as_tibble() %>% 
    mutate(Nstd_ref = Estimate[2],
           Nstd_ref_lower = Q2.5[2],
           Nstd_ref_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Nstd_ref, Nstd_ref_lower, Nstd_ref_upper),
  Sstd_ref_fixef %>% 
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
  Sstd_sens1_fixef %>% 
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
  Sstd_sens3_fixef %>% 
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
  Nstd_sens4_fixef %>% 
    as_tibble() %>% 
    mutate(Nstd_sens4 = Estimate[2],
           Nstd_sens4_lower = Q2.5[2],
           Nstd_sens4_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Nstd_sens4, Nstd_sens4_lower, Nstd_sens4_upper),
  Sstd_sens4_fixef %>% 
    as_tibble() %>% 
    mutate(Sstd_sens4 = Estimate[2],
           Sstd_sens4_lower = Q2.5[2],
           Sstd_sens4_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Sstd_sens4, Sstd_sens4_lower, Sstd_sens4_upper),
  S_PIE_sens4_fixef %>% 
    as_tibble() %>% 
    mutate(S_PIE_sens4 = Estimate[2],
           S_PIE_sens4_lower = Q2.5[2],
           S_PIE_sens4_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(S_PIE_sens4, S_PIE_sens4_lower, S_PIE_sens4_upper),
  # case 11
  Nstd_sens5_fixef %>% 
    as_tibble() %>% 
    mutate(Nstd_sens5 = Estimate[2],
           Nstd_sens5_lower = Q2.5[2],
           Nstd_sens5_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Nstd_sens5, Nstd_sens5_lower, Nstd_sens5_upper),
  Sstd_sens5_fixef %>% 
    as_tibble() %>% 
    mutate(Sstd_sens5 = Estimate[2],
           Sstd_sens5_lower = Q2.5[2],
           Sstd_sens5_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(Sstd_sens5, Sstd_sens5_lower, Sstd_sens5_upper),
  S_PIE_sens5_fixef %>% 
    as_tibble() %>% 
    mutate(S_PIE_sens5 = Estimate[2],
           S_PIE_sens5_lower = Q2.5[2],
           S_PIE_sens5_upper = Q97.5[2]) %>% 
    slice(1) %>% 
    select(S_PIE_sens5, S_PIE_sens5_lower, S_PIE_sens5_upper))
                           

# study level
study_coefs <- left_join(
  # reference fits first
  Nstd_ref_grp_coefs %>% 
    mutate(Nstd_ref = Slope,
           Nstd_ref_upper = Slope_upper, 
           Nstd_ref_lower = Slope_upper) %>% 
    select(dataset_label, Nstd_ref, Nstd_ref_lower, Nstd_ref_upper),
  Sstd_ref_grp_coefs %>% 
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
      Sstd_sens1_grp_coefs %>% 
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
      Sstd_sens3_grp_coefs %>% 
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
      Nstd_sens4_grp_coefs %>% 
        mutate(Nstd_sens4 = Slope,
               Nstd_sens4_upper = Slope_upper, 
               Nstd_sens4_lower = Slope_upper) %>% 
        select(dataset_label, Nstd_sens4, Nstd_sens4_lower, Nstd_sens4_upper),
      by = 'dataset_label') %>% 
    left_join(
      Sstd_sens4_grp_coefs %>% 
        mutate(Sstd_sens4 = Slope,
               Sstd_sens4_upper = Slope_upper, 
               Sstd_sens4_lower = Slope_upper) %>% 
        select(dataset_label, Sstd_sens4, Sstd_sens4_lower, Sstd_sens4_upper),  
      by = 'dataset_label') %>% 
    inner_join(
      S_PIE_sens4_grp_coefs %>% 
        mutate(S_PIE_sens4 = Slope,
               S_PIE_sens4_upper = Slope_upper, 
               S_PIE_sens4_lower = Slope_upper) %>% 
        select(dataset_label, S_PIE_sens4, S_PIE_sens4_lower, S_PIE_sens4_upper),
      by = 'dataset_label') %>% 
    inner_join(
      # case 11
      Nstd_sens5_grp_coefs %>% 
        mutate(Nstd_sens5 = Slope,
               Nstd_sens5_upper = Slope_upper, 
               Nstd_sens5_lower = Slope_upper) %>% 
        select(dataset_label, Nstd_sens5, Nstd_sens5_lower, Nstd_sens5_upper),
      by = 'dataset_label') %>% 
    left_join(
      Sstd_sens5_grp_coefs %>% 
        mutate(Sstd_sens5 = Slope,
               Sstd_sens5_upper = Slope_upper, 
               Sstd_sens5_lower = Slope_upper) %>% 
        select(dataset_label, Sstd_sens5, Sstd_sens5_lower, Sstd_sens5_upper),  
      by = 'dataset_label') %>% 
    inner_join(
      S_PIE_sens5_grp_coefs %>% 
        mutate(S_PIE_sens5 = Slope,
               S_PIE_sens5_upper = Slope_upper, 
               S_PIE_sens5_lower = Slope_upper) %>% 
        select(dataset_label, S_PIE_sens5, S_PIE_sens5_lower, S_PIE_sens5_upper),
      by = 'dataset_label')


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
             aes(x = Sstd_ref, y = Sstd_sens4, colour = 'sens8'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = Sstd_ref, y = Sstd_sens4, colour = 'sens8'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = Sstd_ref_lower, xmax = Sstd_ref_upper, y = Sstd_sens4, colour = 'sens8'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = Sstd_ref, ymin = Sstd_sens4_lower, ymax = Sstd_sens4_upper, colour = 'sens8')
                 ) +
  # reference versus case 11
  geom_point(data = study_coefs,
             aes(x = Sstd_ref, y = Sstd_sens5, colour = 'sens11'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = Sstd_ref, y = Sstd_sens5, colour = 'sens11'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = Sstd_ref_lower, xmax = Sstd_ref_upper, y = Sstd_sens5, colour = 'sens11'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = Sstd_ref, ymin = Sstd_sens5_lower, ymax = Sstd_sens5_upper, colour = 'sens11')) +
  # 1:1 line and zero lines
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  # geom_vline(xintercept = 0, lty = 2) +
  # geom_hline(yintercept = 0, lty = 2) +
  scale_colour_manual(name = 'Imputation',
                      guide = F,
                      values = c('sens1' = '#003f5c',
                                 'sens3' = '#7a5195',
                                 'sens8' = '#ef5675',
                                 'sens11' = '#ffa600'),
                      labels = c('1', '2', '3', '4')) +
  labs(x = '',
       y = '',
       subtitle = expression(paste('Standardised species richness'))) +# (', S[std], ')
  theme_bw() +
  theme(text = element_text(size = 7))


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
             aes(x = Nstd_ref, y = Nstd_sens4, colour = 'sens8'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = Nstd_ref, y = Nstd_sens4, colour = 'sens8'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = Nstd_ref_lower, xmax = Nstd_ref_upper, y = Nstd_sens4, colour = 'sens8'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = Nstd_ref, ymin = Nstd_sens4_lower, ymax = Nstd_sens4_upper, colour = 'sens8')) +
  # reference versus case 11
  geom_point(data = study_coefs,
             aes(x = Nstd_ref, y = Nstd_sens5, colour = 'sens11'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = Nstd_ref, y = Nstd_sens5, colour = 'sens11'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = Nstd_ref_lower, xmax = Nstd_ref_upper, y = Nstd_sens5, colour = 'sens11'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = Nstd_ref, ymin = Nstd_sens5_lower, ymax = Nstd_sens5_upper, colour = 'sens11')) +
  # 1:1 line and zero lines
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  # geom_vline(xintercept = 0, lty = 2) +
  # geom_hline(yintercept = 0, lty = 2) +
  scale_colour_manual(name = 'Imputation',
                      values = c('sens1' = '#003f5c',
                                 'sens3' = '#7a5195',
                                 'sens8' = '#ef5675',
                                 'sens11' = '#ffa600'),
                      labels = c('1', '2', '3', '4')) +
  labs(x = '',
       y = '',
       subtitle = expression(paste('Standardised number of individuals'))) +# (', N[std], ')
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.background = element_blank(),
        # legend.box = element_blank(),
        legend.justification = c(1,0),
        text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6)) +
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
             aes(x = S_PIE_ref, y = S_PIE_sens4, colour = 'sens8'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = S_PIE_ref, y = S_PIE_sens4, colour = 'sens8'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = S_PIE_ref_lower, xmax = S_PIE_ref_upper, y = S_PIE_sens4, colour = 'sens8'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = S_PIE_ref, ymin = S_PIE_sens4_lower, ymax = S_PIE_sens4_upper, colour = 'sens8')) +
  # reference versus case 11
  geom_point(data = study_coefs,
             aes(x = S_PIE_ref, y = S_PIE_sens5, colour = 'sens11'),
             alpha = 0.3) +
  geom_point(data = fixed_effects,
             aes(x = S_PIE_ref, y = S_PIE_sens5, colour = 'sens11'),
             size = 2) +
  geom_errorbarh(data = fixed_effects,
                 aes(xmin = S_PIE_ref_lower, xmax = S_PIE_ref_upper, y = S_PIE_sens5, colour = 'sens11'),
                 height = 0) +
  geom_linerange(data = fixed_effects,
                 aes(x = S_PIE_ref, ymin = S_PIE_sens5_lower, ymax = S_PIE_sens5_upper, colour = 'sens11')) +
  # 1:1 line and zero lines
  scale_colour_manual(name = 'Imputation',
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
       subtitle = expression(paste('Standardised evenness'))) +# (', S[PIE], ')
  theme_bw() +
  theme(text = element_text(size = 7))

# include table showing lack of sensitivity to studies with pooled sampling designs
top1 <- cowplot::ggdraw() + 
  cowplot::draw_image(paste0(path2wd, 'intermediate_results/table_inset.png'),
                      clip = 'on',
                      scale = 0.75) +
  draw_figure_label(label = 'a',
                    size = 8, fontface = 'bold')

bottom <- cowplot::plot_grid(Nstd_robust,
                   Sstd_robust,
                   S_PIE_robust,
                   nrow = 1, align = 'hv',
                   labels = c('b', 'c', 'd'), label_size = 8, label_fontface = 'bold') +
  cowplot::draw_label('Reference imputation',
                      y = 0.03, size = 7) +
  cowplot::draw_label('Alternate imputation',
                      x = 0.01, angle = 90, size = 7)
cowplot::plot_grid(top1, bottom, nrow = 2)
# 2 column
# set local directory
ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision3/figures/Ex_Dat_Fig4.png',
       width = 183,
       height = 120,
       units = 'mm')
