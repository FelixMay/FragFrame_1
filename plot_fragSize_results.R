# code to plot regressions for all metrics (no interaction models)
library(tidyverse)
library(brms)
source('~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/fragSize_coef_wrangle.R')

#---- regression plots showing study-level slopes-----
setwd('~/Dropbox/Frag Database (new)/analysis_apr19/figures/')
# setwd(paste(path2temp,"figs/", sep = ""))

# plot to generate legend 
taxa_legend <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std_2, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = Sstd2_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.75) +
  # scale_colour_viridis_d(guide=F) +
  # scale_color_grey(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  theme_bw() +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        text = element_text(size = 13)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
source('~/Dropbox/1current/R_random/functions/gg_legend.R')
taxa_colour = gg_legend(taxa_legend)

# S_std_2 (Felix recommends this one for the main result)
S_std_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std_2, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = Sstd2_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.75) +
  # fixed effect
  geom_line(data = Sstd2_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1.5) +
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
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(1, 2, 4, 8, 16, 32,64, 128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  # scale_color_grey(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = expression(paste('Species richness (', S[std], ')')),
       tag = 'b') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 13))

# Sn
Sn_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = Sn_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.75) +
  # fixed effect
  geom_line(data = Sn_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1.5) +
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
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(1, 2, 4, 8, 16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  # scale_colour_grey(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = expression(paste('Rarefied richness (', S[n], ')')),
       tag = 'a') +
  theme_bw() +
  theme(legend.position = 'none', 
        text = element_text(size = 13))

# S_PIE
Spie_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_PIE, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = S_PIE_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.75) +
  # fixed effect
  geom_line(data = S_PIE_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1.5) +
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
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(1, 2, 4, 8, 16,32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  # scale_color_grey(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = expression(paste('Evenness (', S[PIE], ')')),
       tag = 'c') +
  theme_bw() +
  theme(legend.position = 'none', text = element_text(size = 13))

# S_cov
ylab = expression(atop('Coverage standardised'), paste('richness (', S[cov], ')'))
Scov_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_cov, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = Scov_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.75) +
  # fixed effect
  geom_line(data = Scov_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1.5) +
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
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(1, 2, 4, 8, 16,32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = expression(atop('Coverage standardised', paste('richness (', S[cov], ')'))),
       tag = 'b') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 13))

# S_cov
Schao_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_chao, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = Schao_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.75) +
  # fixed effect
  geom_line(data = Schao_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1.5) +
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
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(1, 2, 16, 128, 1024)) +#
  # scale_colour_viridis_d(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = expression(paste('Asymptotic richness (', S[chao], ')')),
       tag = 'c') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 13))

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
               size = 0.75) +
  # fixed effect
  geom_line(data = Nstd_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1.5) +
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
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000)) +
  # scale_colour_viridis_d(guide=F) +
  # scale_color_grey(guide = F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = expression(paste('Total abundance (', N[std], ')')),
       tag = 'a') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 13))

# main text: N, S, Spie
top <- cowplot::plot_grid(taxa_colour)
bottom <- cowplot::plot_grid(Nstd_regPlot,
                   S_std_regPlot,
                   Spie_regPlot,
                   nrow = 1, align = 'hv')

cowplot::plot_grid(top, bottom, 
                   nrow = 2,
                   rel_heights = c(0.1,1)) +
  cowplot::draw_label('Fragment size (hectares)', y = 0.05)

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/fig2_taxa_colour.png', width = 250, height = 80, units = 'mm')

bottom_supp <- cowplot::plot_grid(Sn_regPlot,
                                  Scov_regPlot,
                                  Schao_regPlot,
                                  nrow = 1, align = 'hv')
cowplot::plot_grid(top, bottom_supp, 
                   nrow = 2,
                   rel_heights = c(0.1,1)) +
  cowplot::draw_label('Fragment size (hectares)', y = 0.05)
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/figSx_otherMetrics_taxa_color.png', width = 250, height = 80, units = 'mm')

##---coef plots---------
# get the metadata...
meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep = ';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

# create unique dataframes with columns for plotting relationships between change in different metrics
# there are different studies retained for each metric, so put 'em together for each plot separately
S_std_study_slope <- Sstd2_lognorm_fragSize_group_coefs %>% 
    mutate(S_std2_slope = Slope,
           S_std2_lower = Slope_lower,
           S_std2_upper = Slope_upper) %>% 
    select(dataset_label, S_std2_slope, S_std2_lower, S_std2_upper) %>% 
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
             select(dataset_label, S_std2_slope, S_std2_lower, S_std2_upper),
           Nstd_study_slope,
           by = 'dataset_label') %>% 
  # mutate(Sstd_change = ifelse(S_cov_lower > 0, 'delta S_cov > 0', 'no change S_cov '),
  #        N_change = ifelse(S_n_lower > 0, 'delta S_n > 0', 'no change S_n')) %>% 
  ggplot() +
  # facet_wrap(time.since.fragmentation ~ taxa) +
  geom_point(aes(x = N_std_slope, y = S_std2_slope),#, colour = interaction(Scov_change, Sn_change)
             size = 2) +
  geom_linerange(aes(x = N_std_slope, ymin = S_std2_lower, ymax = S_std2_upper),
                 size = 0.25,
                 alpha = 0.2) +
  geom_errorbarh(aes(xmin = N_std_lower, xmax = N_std_upper, y = S_std2_slope ),
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
             select(dataset_label, S_std2_slope, S_std2_lower, S_std2_upper),
           Spie_study_slope,
           by = 'dataset_label') %>% 
  # mutate(Sstd_change = ifelse(S_cov_lower > 0, 'delta S_cov > 0', 'no change S_cov '),
  #        N_change = ifelse(S_n_lower > 0, 'delta S_n > 0', 'no change S_n')) %>% 
  ggplot() +
  # facet_wrap(time.since.fragmentation ~ taxa) +
  geom_point(aes(x = S_PIE_slope, y = S_std2_slope),#, colour = interaction(Scov_change, Sn_change)
             size = 2) +
  geom_linerange(aes(x = S_PIE_slope, ymin = S_std2_lower, ymax = S_std2_upper),
                 size = 0.25,
                 alpha = 0.2) +
  geom_errorbarh(aes(xmin = S_PIE_lower, xmax = S_PIE_upper, y = S_std2_slope ),
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
  facet_wrap(time.since.fragmentation ~ taxa) +
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
  

