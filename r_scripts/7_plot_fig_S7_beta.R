# code to plot beta-diversity regression
# need to execute 0_init_dirs_load_packages.R first

source(paste0(path2wd, 'r_scripts/5d_beta_frag_coef_wrangle.R'))

##-------- plot the regressions----------------
turnover_regression <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Jaccard'),
             aes(x = frag_size_num.x/frag_size_num.y, y = repl, colour = taxa),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Jtu_z1i_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin.x/xmin.y,
                   # trick: max x-value requires a min in the denominator 
                   xend = xmax.x/xmin.y,
                   y = plogis(Intercept + Slope * cxmin),
                   yend = plogis(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Jtu_z1i_fitted, 
            aes(x = frag_size_num.x/frag_size_num.y,
                y = plogis(Estimate)),
            size = 1) +
  # fixed effect uncertainty
  geom_ribbon(data = Jtu_z1i_fitted,
              aes(x = frag_size_num.x/frag_size_num.y,
                  ymin = plogis(Q2.5),
                  ymax = plogis(Q97.5)),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = 'Turnover (Jaccard)',
       tag = 'a') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 7),
        plot.tag = element_text(size = 8, face = 'bold'))


turnover_regression_ruzicka <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Ruzicka'),
             aes(x = frag_size_num.x/frag_size_num.y, y = repl, colour = taxa),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Rtu_z1i_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin.x/xmin.y,
                   # trick: max x-value requires a min in the denominator 
                   xend = xmax.x/xmin.y,
                   y = plogis(Intercept + Slope * cxmin),
                   yend = plogis(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect for Ruzicka
  geom_line(data = Rtu_z1i_fitted,
            aes(x = frag_size_num.x/frag_size_num.y,
                y = plogis(Estimate)),
            size = 1) +
  # fixed effect uncertainty: Ruzicka
  geom_ribbon(data = Rtu_z1i_fitted,
              aes(x = frag_size_num.x/frag_size_num.y,
                  ymin = plogis(Q2.5),
                  ymax = plogis(Q97.5)),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = 'Balanced abundance (Ruzicka)',
       tag = 'b') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 7),
        plot.tag = element_text(size = 8, face = 'bold'))

nestedness_regression <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Jaccard'),
             aes(x = frag_size_num.x/frag_size_num.y, y = rich, colour = taxa),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Jne_zi_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin.x/xmin.y,
                   xend = xmax.x/xmin.y,
                   y = plogis(Intercept + Slope * cxmin),
                   yend = plogis(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Jne_zi_fitted, 
            aes(x = frag_size_num.x/frag_size_num.y,
                y = plogis(Estimate)),
            size = 1) +
  # fixed effect uncertainty
  geom_ribbon(data = Jne_zi_fitted,
              aes(x = frag_size_num.x/frag_size_num.y,
                  ymin = plogis(Q2.5),
                  ymax = plogis(Q97.5)),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = 'Nestedness (Jaccard)', 
       tag = 'c') +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 7),
        plot.tag = element_text(size = 8, face = 'bold'))

nestedness_regression_ruzicka <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Ruzicka'),
             aes(x = frag_size_num.x/frag_size_num.y, y = rich, colour = taxa),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Rne_zi_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin.x/xmin.y,
                   xend = xmax.x/xmin.y,
                   y = plogis(Intercept + Slope * cxmin),
                   yend = plogis(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Rne_zi_fitted,
            aes(x = frag_size_num.x/frag_size_num.y,
                y = plogis(Estimate)),
            size = 1) +
  # fixed effect uncertainty: Ruzicka
  geom_ribbon(data = Rne_zi_fitted,
              aes(x = frag_size_num.x/frag_size_num.y,
                  ymin = plogis(Q2.5),
                  ymax = plogis(Q97.5)),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = 'Abundance gradient (Ruzicka)',
       tag = 'd') +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 7),
        plot.tag = element_text(size = 8, face = 'bold'))

# plot to generate legend 
taxa_legend <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Ruzicka'),
             aes(x = frag_size_num.x/frag_size_num.y, y = rich, colour = taxa),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Rne_zi_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin.x/xmin.y,
                   xend = xmax.x/xmin.y,
                   y = plogis(Intercept + Slope * cxmin),
                   yend = plogis(Intercept + Slope * cxmax)),
               size = 0.5) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  theme_bw() +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        text = element_text(size = 7)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1)))

source(paste0(path2wd, 'r_scripts/99_gg_legend.R'))
taxa_colour = gg_legend(taxa_legend)

top <- cowplot::plot_grid(taxa_colour)
bottom <- cowplot::plot_grid(turnover_regression,
                             turnover_regression_ruzicka,
                             nestedness_regression,
                             nestedness_regression_ruzicka,
                   nrow = 2, align = 'hv')

cowplot::plot_grid(top, bottom, 
                   nrow = 2,
                   rel_heights = c(0.075,1)) +
  cowplot::draw_label('Ratio of fragment sizes (log-scale)',
                      y = 0.02, size = 7)

# plot for 2 column width
ggsave(paste0(path2wd, 'extended_data_figs_tabs/Ex_Dat_Fig7.png'),
       width = 183,
       height = 170,
       units = 'mm')
