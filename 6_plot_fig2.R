# need to execute 0_init_dirs_load_packages.R first

# code to plot regressions for all metrics (no interaction models)

# get the coefficients for all the results
source(paste0(path2wd, '5a_fragSize_coef_wrangle.R'))

#---- regression plots showing study-level slopes-----
setwd(paste0(path2Dropbox, '/Manuscript for Nature/revision3/figures/'))
# setwd(paste(path2temp,"figs/", sep = ""))

# plot to generate legend 
taxa_legend <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std2_mean, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = Sstd2_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # scale_colour_viridis_d(guide=F) +
  # scale_color_grey(guide=F) +
  scale_colour_brewer(name = 'Taxon group', type = 'qual', palette = 'Dark2',
                      labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))+
  theme_bw() +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        text = element_text(size = 8)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1)))

source('~/Dropbox/1current/R_random/functions/gg_legend.R')
taxa_colour = gg_legend(taxa_legend)

# S_std (Felix recommends this one for the main result)
S_std_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std2_mean, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = Sstd2_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Sstd2_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1) +
  # fixed effect uncertainty
  geom_ribbon(data = Sstd2_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.01, y = Inf, hjust = 0.1, vjust = 1.4,
           label = paste("beta == ", #[Frag.~size]
                         round(Sstd2_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         " (",
                         round(Sstd2_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
                         " - ",
                         round(Sstd2_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
                         ")"),  
           parse = T, size = 2) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(1, 2, 4, 8, 16, 32,64, 128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  # scale_color_grey(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = expression(paste('Standardised species richness')),
       tag = 'c') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'))

# Sn
Sn_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n_mean, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = Sn_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Sn_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1) +
  # fixed effect uncertainty
  geom_ribbon(data = Sn_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.01, y = Inf, hjust = 0.1, vjust = 1.4,
           label = paste("beta == ", #[Frag.~size]
                         round(Sn_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         "  (",
                         round(Sn_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
                         " - ",
                         round(Sn_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
                         ")"),  
           parse = T, size = 2) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(1, 2, 4, 8, 16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  # scale_colour_grey(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = expression(paste('Rarefied richness')),
       tag = 'a') +
  theme_bw() +
  theme(legend.position = 'none', 
        text = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'))

# S_PIE
Spie_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_PIE_mean, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = S_PIE_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = S_PIE_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1) +
  # fixed effect uncertainty
  geom_ribbon(data = S_PIE_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.01, y = Inf, hjust = 0.1, vjust = 1.4,
           label = paste("beta == ", #[Frag.~size]
                         round(S_PIE_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         "  (",
                         round(S_PIE_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
                         " - ",
                         round(S_PIE_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
                         ")"),  
           parse = T, size = 2) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(1, 2, 4, 8, 16,32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  # scale_color_grey(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = expression(paste('Standardised evenness')),
       tag = 'd') +
  theme_bw() +
  theme(legend.position = 'none', 
        text = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'))

# S_cov
Scov_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_cov_mean, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = Scov_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Scov_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1) +
  # fixed effect uncertainty
  geom_ribbon(data = Scov_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.01, y = Inf, hjust = 0.1, vjust = 1.4,
           label = paste("beta == ", #[Frag.~size]
                         round(Scov_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         "  (",
                         round(Scov_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
                         " - ",
                         round(Scov_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
                         ")"),  
           parse = T, size = 2) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(1, 2, 4, 8, 16,32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = 'Coverage standardised richness',
       tag = 'b') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'))

# S_cov
Schao_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_chao_mean, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = Schao_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Schao_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1) +
  # fixed effect uncertainty
  geom_ribbon(data = Schao_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.01, y = Inf, hjust = 0.1, vjust = 1.4,
           label = paste("beta == ", #[Frag.~size]
                         round(Schao_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         "  (",
                         round(Schao_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
                         " - ",
                         round(Schao_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
                         ")"),  
           parse = T, size = 2) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(1, 2, 16, 128, 1024)) +#
  # scale_colour_viridis_d(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = expression(paste('Asymptotic richness')),
       tag = 'c') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'))

# N_std
Nstd_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = N_std, 
                 colour = taxa
                 ),
             size = 1, alpha = 0.25) +
  geom_segment(data = Nstd_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Nstd_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1) +
  # fixed effect uncertainty
  geom_ribbon(data = Nstd_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.01, y = Inf, hjust = 0.1, vjust = 1.4,
           label = paste("beta == ", #[Frag.~size]
                         round(Nstd_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         "  (",
                         round(Nstd_lognorm_fragSize_fixef['c.lfs','Q2.5'],3),
                         " - ",
                         round(Nstd_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
                         ")"),  
           parse = T, size = 2) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_colour_viridis_d(guide=F) +
  # scale_color_grey(guide = F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = expression(paste('Standardised number of individuals')),
       tag = 'b') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'))

# main text: N, S, Spie
top <- cowplot::plot_grid(taxa_colour)
bottom <- cowplot::plot_grid(Nstd_regPlot,
                   S_std_regPlot,
                   Spie_regPlot,
                   nrow = 1, align = 'hv') +
  cowplot::draw_label('Fragment size (hectares)', y = 0.05, size = 9)

# bottom_panel <- cowplot::plot_grid(top, bottom, 
#                    nrow = 2,
#                    rel_heights = c(0.1,1)) 

# new figure two has map on it
source(paste0(path2wd, '6_plot_map.R'))

cowplot::plot_grid(map_taxa,
                   bottom,
                   nrow = 2)
# plots wrangled for two column print 
ggsave('fig2_2column.pdf', width = 183, height = 120, units = 'mm')

bottom_supp <- cowplot::plot_grid(Sn_regPlot,
                                  Scov_regPlot,
                                  Schao_regPlot,
                                  nrow = 1, align = 'hv')
cowplot::plot_grid(top, bottom_supp, 
                   nrow = 2,
                   rel_heights = c(0.1,1)) +
  cowplot::draw_label('Fragment size (hectares)', y = 0.05)
# ggsave('figS1_otherMetrics_revision.png', width = 250, height = 80, units = 'mm')

##---coef plots---------
# get the metadata...
meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep = ';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

meta <- meta %>% 
  separate(coordinates, into = c('y', 'x'), sep = ', ', remove = F) %>% 
  mutate(x = as.numeric(x),
         y = as.numeric(y),
         latitude = climate)

# create unique dataframes with columns for plotting relationships between change in different metrics
# there are different studies retained for each metric, so put 'em together for each plot separately
S_std_study_slope <- Sstd2_lognorm_fragSize_group_coefs %>% 
    mutate(S_std_slope = Slope,
           error = error,
           S_std_lower = Slope_lower,
           S_std_upper = Slope_upper) %>% 
    select(dataset_label, S_std_slope, error, S_std_lower, S_std_upper) %>% 
  left_join(meta,
            by = 'dataset_label')

Sn_study_slope<- Sn_lognorm_fragSize_group_coefs %>% 
    mutate(S_n_slope = Slope,
           S_n_lower = Slope_lower,
           S_n_upper = Slope_upper) %>% 
    select(dataset_label, S_n_slope, S_n_lower, S_n_upper) %>% 
  left_join(meta,
            by = 'dataset_label')

Scov_study_slope<- Scov_lognorm_fragSize_group_coefs %>% 
  mutate(S_cov_slope = Slope,
         S_cov_lower = Slope_lower,
         S_cov_upper = Slope_upper) %>% 
  select(dataset_label, S_cov_slope, S_cov_lower, S_cov_upper) %>% 
  left_join(meta,
            by = 'dataset_label')

Schao_study_slope<- Schao_lognorm_fragSize_group_coefs %>% 
  mutate(S_chao_slope = Slope,
         S_chao_lower = Slope_lower,
         S_chao_upper = Slope_upper) %>% 
  select(dataset_label, S_chao_slope, S_chao_lower, S_chao_upper) %>% 
  left_join(meta,
            by = 'dataset_label')

Spie_study_slope<- S_PIE_fragSize_group_coefs %>% 
  mutate(S_PIE_slope = Slope,
         S_PIE_lower = Slope_lower,
         S_PIE_upper = Slope_upper) %>% 
  select(dataset_label, S_PIE_slope, S_PIE_lower, S_PIE_upper) %>% 
  left_join(meta,
            by = 'dataset_label')

Nstd_study_slope <- Nstd_fragSize_group_coefs %>% 
  mutate(N_std_slope = Slope,
         N_std_lower = Slope_lower,
         N_std_upper = Slope_upper) %>% 
  select(dataset_label, N_std_slope, N_std_lower, N_std_upper) %>% 
  left_join(meta,
            by = 'dataset_label')

# delta S_std ~ delta N
inner_join(S_std_study_slope %>% 
             select(dataset_label, S_std_slope, S_std_lower, S_std_upper),
           Nstd_study_slope,
           by = 'dataset_label') %>% 
  # mutate(Sstd2_change = ifelse(S_cov_lower > 0, 'delta S_cov > 0', 'no change S_cov '),
  #        N_change = ifelse(S_n_lower > 0, 'delta S_n > 0', 'no change S_n')) %>% 
  ggplot() +
  # facet_wrap(time.since.fragmentation ~ taxa) +
  geom_point(aes(x = N_std_slope, y = S_std_slope),#, colour = interaction(Scov_change, Sn_change)
             size = 2) +
  geom_linerange(aes(x = N_std_slope, ymin = S_std_lower, ymax = S_std_upper),
                 size = 0.25,
                 alpha = 0.2) +
  geom_errorbarh(aes(xmin = N_std_lower, xmax = N_std_upper, y = S_std_slope ),
                 size = 0.25,
                 alpha = 0.2) +
  geom_abline(intercept = 0, 
              slope = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

# delta S_std ~ delta S_PIE
inner_join(S_std_study_slope %>% 
             select(dataset_label, S_std_slope, S_std_lower, S_std_upper),
           Spie_study_slope,
           by = 'dataset_label') %>% 
  # mutate(Sstd2_change = ifelse(S_cov_lower > 0, 'delta S_cov > 0', 'no change S_cov '),
  #        N_change = ifelse(S_n_lower > 0, 'delta S_n > 0', 'no change S_n')) %>% 
  ggplot() +
  # facet_wrap(time.since.fragmentation ~ taxa) +
  geom_point(aes(x = S_PIE_slope, y = S_std_slope),#, colour = interaction(Scov_change, Sn_change)
             size = 2) +
  geom_linerange(aes(x = S_PIE_slope, ymin = S_std_lower, ymax = S_std_upper),
                 size = 0.25,
                 alpha = 0.2) +
  geom_errorbarh(aes(xmin = S_PIE_lower, xmax = S_PIE_upper, y = S_std_slope ),
                 size = 0.25,
                 alpha = 0.2) +
  geom_abline(intercept = 0, 
              slope = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

# what about Scov ~ Sn: difference between these two is due to deltaN
inner_join(Sn_study_slope %>% 
             select(dataset_label, S_n_slope, S_n_lower, S_n_upper),
           Scov_study_slope,
           by = 'dataset_label') %>% 
  mutate(Scov_change = ifelse(S_cov_lower > 0, 'delta S_cov > 0', 'no change S_cov '),
         Sn_change = ifelse(S_n_lower > 0, 'delta S_n > 0', 'no change S_n')) %>% 
  ggplot() +
  # facet_wrap(time.since.fragmentation ~ taxa) +
  geom_point(aes(x = S_n_slope, y = S_cov_slope, colour = interaction(Scov_change, Sn_change)),
             size = 2) +
  geom_linerange(aes(x = S_n_slope, ymin = S_cov_lower, ymax = S_cov_upper),
                 size = 0.25,
                 alpha = 0.2) +
  geom_errorbarh(aes(xmin = S_n_lower, xmax = S_n_upper, y = S_cov_slope ),
                 size = 0.25,
                 alpha = 0.2) +
  geom_abline(intercept = 0, 
              slope = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
  

ggplot() +
  geom_point(data = S_std_study_slope,
             aes(x = abs(y), y = S_std_slope),
             size = 1.5) +
  geom_linerange(data = S_std_study_slope,
                 aes(x = abs(y), ymin = S_std_lower, ymax = S_std_upper),
                 size = 0.5, alpha = 0.5) +
  stat_smooth(data = S_std_study_slope,
              aes(x = abs(y), y = S_std_slope),
              method = 'lm')

# ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision1/figures/decay_latitude.png', width = 120, height = 120, units = 'mm')
with(S_std_study_slope, summary(lm(S_std_slope ~ abs(y))))
# include uncertainty in the study-level estimate of the 
S_std_study_slope <- S_std_study_slope %>% 
  mutate(abs_lat = abs(y))

slope_lat <- brms::brm(bf(S_std_slope | se(error) ~ abs_lat + (1|dataset_label)),
                       data = S_std_study_slope,
                       cores = 4, chains = 4)

pp_check(slope_lat)
plot(slope_lat)

slope_lat_fixef <- cbind(slope_lat$data,
                             fitted(slope_lat, 
                                    re_formula = NA)) %>% 
  as_tibble()

ggplot() +
  geom_point(data = S_std_study_slope,
             aes(x = abs(y), y = S_std_slope),
             size = 1.5) +
  geom_linerange(data = S_std_study_slope,
                 aes(x = abs(y), ymin = S_std_lower, ymax = S_std_upper),
                 size = 0.5, alpha = 0.5) +
  geom_line(data = slope_lat_fixef,
            aes(x = abs_lat, y = Estimate), 
            size = 1.2) +
  geom_ribbon(data = slope_lat_fixef,
              aes(x = abs_lat, ymin = Q2.5, ymax = Q97.5), 
              alpha = 0.5) +
  labs(x = '|Latitude|',
       y = 'Richness ~ fragment size slope estimate')
              
# ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision1/figures/decay_latitude_option1.png',
#        width = 120, height = 120, units = 'mm')              

slope_lat_taxa <- brms::brm(bf(S_std_slope | se(error) ~ abs_lat*taxa + (1|dataset_label)),
                       data = S_std_study_slope,
                       cores = 4, chains = 4)

slope_lat <- add_criterion(slope_lat, criterion = c('loo', 'waic'))
slope_lat_taxa <- add_criterion(slope_lat_taxa, criterion = c('loo', 'waic'))

model_weights(slope_lat, slope_lat_taxa)
loo_compare(slope_lat, slope_lat_taxa, criterion = 'waic')
