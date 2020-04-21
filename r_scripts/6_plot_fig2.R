# need to execute 0_init_dirs_load_packages.R first

# code to plot regressions for all metrics (no interaction models)

# get the coefficients for all the results
source(paste0(path2wd, 'r_scripts/5a_fragSize_coef_wrangle.R'))

#---- regression plots showing study-level slopes----

# plot to generate legend 
taxa_legend <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std_mean, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = Sstd_lognorm_fragSize_group_coefs,
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
        text = element_text(size = 7)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1)))

# function to plot legend only
source(paste0(path2wd, 'r_scripts/99_gg_legend.R'))
taxa_colour = gg_legend(taxa_legend)

# S_std: standardised species richness
S_std_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std_mean, colour = taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = Sstd_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Sstd_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1) +
  # fixed effect uncertainty
  geom_ribbon(data = Sstd_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.01, y = Inf, hjust = 0.1, vjust = 1.4,
           label = paste("beta == ", #[Frag.~size]
                         round(Sstd_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         " (",
                         round(Sstd_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
                         " - ",
                         round(Sstd_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
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
        text = element_text(size = 7),
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
        text = element_text(size = 7),
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
        text = element_text(size = 7),
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
        text = element_text(size = 7),
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
        text = element_text(size = 7),
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
        text = element_text(size = 7),
        plot.tag = element_text(size = 8, face = 'bold'))

# main text: N, S, Spie
top <- cowplot::plot_grid(taxa_colour)
bottom <- cowplot::plot_grid(Nstd_regPlot,
                   S_std_regPlot,
                   Spie_regPlot,
                   nrow = 1, align = 'hv') +
  cowplot::draw_label('Fragment size (hectares)', y = 0.05, size = 7)

# run script to plot map
source(paste0(path2wd, 'r_scripts/6_plot_fig2_map.R'))

cowplot::plot_grid(map_taxa,
                   bottom,
                   nrow = 2)

# plots sized for two column print
# set local_directory
# ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision3/figures/fig2_2column.pdf', width = 183, height = 120, units = 'mm')

bottom_supp <- cowplot::plot_grid(Sn_regPlot,
                                  Scov_regPlot,
                                  Schao_regPlot,
                                  nrow = 1, align = 'hv')
cowplot::plot_grid(top, bottom_supp, 
                   nrow = 2,
                   rel_heights = c(0.1,1)) +
  cowplot::draw_label('Fragment size (hectares)', y = 0.05, size = 7)

# set local directory and save
# ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision3/figures/Ex_Dat_Fig2.png', 
#        width = 183, height = 60, units = 'mm')
# 
