# need to execute 0_init_dirs_load_packages.R first

# COLOUR VERSION

# get study-level posterior samples from model with only study-level slope variation
# code to plot each of the groups we are interested in together on one figure: 
# taxa, continent, time, matrix

source(paste0(path2wd, 'r_scripts/5aa_alpha_frag_posterior_wrangle.R'))

## dummy plot for creating separate legend
three_grey_legend <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = Sstd_posterior,
                               aes(x = S_std + Sstd_global,
                                   y = time.since.fragmentation,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9,
                               linetype = 0) +
  theme_bw() +
  labs(y = 'Time since fragmentation',
       x = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363'),
                    labels = c('< 5%', '< 45%',  '50%')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        # legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.text = element_text(size = 5, face = 'plain'),
        legend.title = element_text(size = 6, face = 'plain'),
        legend.margin = margin(),
        legend.box.spacing = unit(c(0,0,0,0), units = 'mm'),
        legend.key.size = unit(2, units = 'mm'),
        plot.margin = unit(c(0,0,0,0), units = 'mm')) #+

source(paste0(path2wd, 'r_scripts/99_gg_legend.R'))
legend <- gg_legend(three_grey_legend)

Sstd_study_posterior_time <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = Sstd_posterior,
                               aes(x = S_std + Sstd_global,
                                   y = time.since.fragmentation,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Sstd_global, probs = 0.025),
                        upper = quantile(Sstd_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.5) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Sstd_global)),
             size = 0.5,
             alpha = 0.8) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5) +
  geom_point(data = Sstd_posterior,
             aes(x = S_std + Sstd_global, 
                 y = time.since.fragmentation),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.25, shape = 18, colour = 'black') +
  geom_text(data = Sstd_posterior %>%
              group_by(time.since.fragmentation) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(time.since.fragmentation, n_study, .keep_all = T),
            aes(x=0.3, y=time.since.fragmentation, 
                label=paste('n[study] == ', n_study)),
            size=1.5,
            nudge_y = 0.15, parse = T) +
  # geom_text(data = Sstd_posterior %>%
  #             group_by(time.since.fragmentation) %>% 
  #             summarise(min_x = min(S_std)),
  #           aes(x=-0.3, y=time.since.fragmentation, 
  #               label=time.since.fragmentation),
  #           size=2.5,
  #           nudge_y = 0.1, 
  #           #parse = T
  #           ) +
  coord_equal(ratio = 0.8/3) +
  theme_bw() +
  labs(y = 'Time since fragmentation',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'c',
       subtitle = 'Time since fragmentation'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                                          '#969696', '#cccccc')) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(vjust = 0, hjust = 1),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 7, hjust = -1.2))#+


Sstd_study_posterior_matrix <- ggplot() +
  geom_density_ridges_gradient(data = Sstd_posterior,
                               aes(x = S_std + Sstd_global,
                                   y = Matrix.category,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Sstd_global, probs = 0.025),
                        upper = quantile(Sstd_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.5) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Sstd_global)),
             size = 0.5,
             alpha = 0.8) +
  geom_vline(xintercept = 0, lty = 2,
             size = 0.5) +
  geom_point(data = Sstd_posterior,
             aes(x = S_std + Sstd_global, 
                 y = Matrix.category),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.25, shape = 18, colour = 'black') +
  geom_text(data = Sstd_posterior %>%
              group_by(Matrix.category) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(Matrix.category, n_study, .keep_all = T),
            aes(x=0.3, y=Matrix.category, 
                label=paste('n[study] == ', n_study)),
            size=1.5,
            nudge_y = 0.15, parse = T) +
  # geom_text(data = Sstd_posterior %>%
  #             group_by(Matrix.category) %>% 
  #             summarise(min_x = min(S_std)),
  #           aes(x=-0.3, y=Matrix.category, 
  #               label=Matrix.category),
  #           size=2.5,
  #           nudge_y = 0.1, 
  #           #parse = T
  # ) +
  theme_bw() +
  coord_equal(ratio = 0.8/3) +
  labs(y = 'Matrix filter',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'd',
       subtitle = 'Matrix filter'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                               '#969696', '#cccccc')) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(vjust = 0, hjust = 1),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 7, hjust = -0.5))

Sstd_study_posterior_taxa <- ggplot() +
  geom_density_ridges_gradient(data = Sstd_posterior,
                               aes(x = S_std + Sstd_global,
                                   y = taxa,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Sstd_global, probs = 0.025),
                        upper = quantile(Sstd_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.5) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Sstd_global)),
             size = 0.5,
             alpha = 0.8) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5) +
  geom_point(data = Sstd_posterior,
             aes(x = S_std + Sstd_global, 
                 y = taxa),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.25, shape = 18, colour = 'black') +
  geom_text(data = Sstd_posterior %>%
              group_by(taxa) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(taxa, n_study, .keep_all = T),
            aes(x=0.3, y=taxa, 
                label=paste('n[study] == ', n_study)),
            size=1.5,
            nudge_y = 0.15, parse = T) +
  # geom_text(data = Sstd_posterior %>%
  #             group_by(taxa) %>% 
  #             summarise(min_x = min(S_std)),
  #           aes(x=-0.3, y=taxa, 
  #               label=taxa),
  #           size=2,
  #           nudge_y = 0.1, 
  #           #parse = T
  # ) +
  theme_bw() +
  coord_equal(ratio = 0.8/5) +
  labs(y = 'Taxon group',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'),
       tag = 'a',
       subtitle = 'Taxon group'
       
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                               '#969696', '#cccccc')) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(vjust = 0, hjust = 1),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 7, hjust = -0.52))

Sstd_study_posterior_continent <- ggplot() +
  geom_density_ridges_gradient(data = Sstd_posterior,
                               aes(x = S_std + Sstd_global,
                                   y = continent8,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Sstd_global, probs = 0.025),
                        upper = quantile(Sstd_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.5) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Sstd_global)),
             size = 0.5,
             alpha= 0.8) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5) +
  geom_point(data = Sstd_posterior,
             aes(x = S_std + Sstd_global, 
                 y = continent8),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.25, shape = 18, colour = 'black') +
  geom_text(data = Sstd_posterior %>%
              group_by(continent8) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(continent8, n_study, .keep_all = T),
            aes(x=0.3, y=continent8, 
                label=paste('n[study] == ', n_study)),
            size=1.5,
            nudge_y = 0.15, parse = T) +
  theme_bw() +
  coord_equal(ratio = 0.8/7) +
  labs(y = 'Continent',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'b',
       subtitle = 'Continent'
  ) +
  scale_y_discrete(
    labels = scales::wrap_format(12),
    expand = c(0.05,0,0.1,0)
  ) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                               '#969696', '#cccccc')) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        # plot.tag.position = 'top',
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(vjust = 0, hjust = 1),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 7, hjust = -0.5))


top = cowplot::plot_grid(NULL, legend, NULL,
                         nrow = 1, rel_widths = c(0.2,0.6,0.2))
# top = cowplot::plot_grid(legend)
# top_column = cowplot::plot_grid(legend)
bottom = cowplot::plot_grid(Sstd_study_posterior_taxa,
                            Sstd_study_posterior_continent,
                            Sstd_study_posterior_time, 
                            Sstd_study_posterior_matrix, 
                            nrow = 2,
                            align = 'hv')
cowplot::plot_grid(top, 
                   bottom, 
                   # top_column,
                   # rel_widths = c(1,0.2),
                   rel_heights = c(0.025,1),
                   nrow = 2) +
  cowplot::draw_label(expression(paste('Standardised species richness ~ fragment size slope estimate')), y = 0.01, size = 7)

# set local directory
# plot sized for for 2 column width 
ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision3/figures/test3_120wide_grey.pdf',
       width = 120,
       height = 110,
       units = 'mm')
# 
##repeat for Nstd and Nstd for supplement
N_continent <- ggplot() +
  geom_density_ridges_gradient(data = Nstd_posterior,
                               aes(x = Nstd + Nstd_global_slope,
                                   y = continent8,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Nstd_global, probs = 0.025),
                        upper = quantile(Nstd_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Nstd_global)),
             size = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5) +
  geom_point(data = Nstd_posterior,
             aes(x = Nstd + Nstd_global_slope, 
                 y = continent8),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.25, shape = 18) +
  geom_text(data = Nstd_posterior %>%
              group_by(continent8) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(continent8, n_study, .keep_all = T),
            aes(x=-0.3, y=continent8, 
                label=paste('n[study] == ', n_study)),
            size=1.5,
            nudge_y = 0.15, parse = T) +
  theme_bw() +
  labs(y = 'Continent',
       x = '',#,#expression(paste('Study-level slope')),
       subtitle = 'Continent',
       tag = 'b'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                               '#969696', '#cccccc')) +
  coord_equal(ratio = 0.8/7) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(vjust = 0, hjust = 1),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 7, hjust = -0.26)) #+

N_taxa <- ggplot() +
  geom_density_ridges_gradient(data = Nstd_posterior,
                               aes(x = Nstd + Nstd_global_slope,
                                   y = taxa,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Nstd_global, probs = 0.025),
                        upper = quantile(Nstd_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Nstd_global)),
             size = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5) +
  geom_point(data = Nstd_posterior,
             aes(x = Nstd + Nstd_global_slope, 
                 y = taxa),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.25, shape = 18) +
  geom_text(data = Nstd_posterior %>%
              group_by(taxa) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(taxa, n_study, .keep_all = T),
            aes(x=-0.3, y=taxa, 
                label=paste('n[study] == ', n_study)),
            size=1.5,
            nudge_y = 0.15, parse = T) +
  theme_bw() +
  labs(y = 'Taxon group',
       x = '',#,#expression(paste('Study-level slope')),
       subtitle = 'Taxon group',
       tag = 'a'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                               '#969696', '#cccccc')) +
  coord_equal(ratio = 0.8/5) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(vjust = 0, hjust = 1),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 7, hjust = -0.5)) #+

N_time <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = Nstd_posterior,
                               aes(x = Nstd + Nstd_global_slope,
                                   y = time.since.fragmentation,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Nstd_global, probs = 0.025),
                        upper = quantile(Nstd_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Nstd_global)),
             size = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5) +
  geom_point(data = Nstd_posterior,
             aes(x = Nstd + Nstd_global_slope, 
                 y = time.since.fragmentation),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.25, shape = 18) +
  geom_text(data = Nstd_posterior %>%
              group_by(time.since.fragmentation) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(time.since.fragmentation, n_study, .keep_all = T),
            aes(x=-0.3, y=time.since.fragmentation, 
                label=paste('n[study] == ', n_study)),
            size=1.5,
            nudge_y = 0.15, parse = T) +
  theme_bw() +
  labs(y = 'Time since fragmentation',
       x = '',#,#expression(paste('Study-level slope')),
       subtitle = 'Time since fragmentation',
       tag = 'c'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                               '#969696', '#cccccc')) +
  coord_equal(ratio = 0.8/3) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(vjust = 0, hjust = 1),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 7, hjust = -0.65)) #+


N_matrix <- ggplot() +
  geom_density_ridges_gradient(data = Nstd_posterior,
                               aes(x = Nstd + Nstd_global_slope,
                                   y = Matrix.category,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Nstd_global, probs = 0.025),
                        upper = quantile(Nstd_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Nstd_global)),
             size = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5) +
  geom_point(data = Nstd_posterior,
             aes(x = Nstd + Nstd_global_slope, 
                 y = Matrix.category),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.25, shape = 18) +
  geom_text(data = Nstd_posterior %>%
              group_by(Matrix.category) %>%
              mutate(n_study = n_distinct(dataset_label)) %>%
              ungroup() %>%
              distinct(Matrix.category, n_study, .keep_all = T),
            aes(x=-0.3, y=Matrix.category,
                label=paste('n[study] == ', n_study)),
            size=1.5,
            nudge_y = 0.15, parse = T) +
  theme_bw() +
  labs(y = 'Matrix filter',
       x = '',#,#expression(paste('Study-level slope')),
       subtitle = 'Matrix filter',
       tag = 'd'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                               '#969696', '#cccccc')) +
  coord_equal(ratio = 0.8/3) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(vjust = 0, hjust = 1),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 7, hjust = -0.42)) #+

bottom1 = cowplot::plot_grid(N_taxa,
                             N_continent,
                             N_time,  
                             N_matrix,
                             nrow = 2
) +
  cowplot::draw_label(expression(paste('Standardised number of individuals ~ fragment size slope estimate')), 
                      y = 0.015, size = 7)

##repeat for S_PIE for supplement
continent <- ggplot() +
  geom_density_ridges_gradient(data = Spie_posterior,
                               aes(x = S_PIE + S_PIE_global,
                                   y = continent8,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(S_PIE_global, probs = 0.025),
                        upper = quantile(S_PIE_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_vline(data = frag_global,
             aes(xintercept = median(S_PIE_global)), 
             size = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5) +
  geom_point(data = Spie_posterior,
             aes(x = S_PIE + S_PIE_global, 
                 y = continent8),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.5, shape = 18) +
  geom_text(data = Spie_posterior %>%
              group_by(continent8) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(continent8, n_study, .keep_all = T),
            aes(x=-0.175, y=continent8, 
                label=paste('n[study] == ', n_study)),
            size=1.5,
            nudge_y = 0.15, parse = T) +
  theme_bw() +
  labs(y = 'Continent',
       x = '',#,#expression(paste('Study-level slope')),
       subtitle = 'Continent',
       tag = 'f'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                               '#969696', '#cccccc')) +
  coord_equal(ratio = 0.5/7) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(vjust = 0, hjust = 1),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 7, hjust = -0.25)) #+

taxa <- ggplot() +
  geom_density_ridges_gradient(data = Spie_posterior,
                               aes(x = S_PIE + S_PIE_global,
                                   y = taxa,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(S_PIE_global, probs = 0.025),
                        upper = quantile(S_PIE_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_vline(data = frag_global,
             aes(xintercept = median(S_PIE_global)),
             size = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5) +
  geom_point(data = Spie_posterior,
             aes(x = S_PIE + S_PIE_global, 
                 y = taxa),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.25, shape = 18) +
  geom_text(data = Spie_posterior %>%
              group_by(taxa) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(taxa, n_study, .keep_all = T),
            aes(x=-0.175, y=taxa, 
                label=paste('n[study] == ', n_study)),
            size=1.5,
            nudge_y = 0.15, parse = T) +
  theme_bw() +
  labs(y = 'Taxon group',
       x = '',#,#expression(paste('Study-level slope')),
       subtitle = 'Taxon group',
       tag = 'e'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                               '#969696', '#cccccc')) +
  coord_equal(ratio = 0.5/5) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(vjust = 0, hjust = 1),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 7, hjust = -0.5)) #+

time <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = Spie_posterior,
                               aes(x = S_PIE + S_PIE_global,
                                   y = time.since.fragmentation,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(S_PIE_global, probs = 0.025),
                        upper = quantile(S_PIE_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_vline(data = frag_global,
             aes(xintercept = median(S_PIE_global)),
             size = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5) +
  geom_point(data = Spie_posterior,
             aes(x = S_PIE + S_PIE_global, 
                 y = time.since.fragmentation),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.25, shape = 18) +
  geom_text(data = Spie_posterior %>%
              group_by(time.since.fragmentation) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(time.since.fragmentation, n_study, .keep_all = T),
            aes(x=-0.175, y=time.since.fragmentation, 
                label=paste('n[study] == ', n_study)),
            size=1.5,
            nudge_y = 0.15, parse = T) +
  theme_bw() +
  labs(y = 'Time since fragmentation',
       x = '',#,#expression(paste('Study-level slope')),
       subtitle = 'Time since fragmentation',
       tag = 'g'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                               '#969696', '#cccccc')) +
  coord_equal(ratio = 0.5/3) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(vjust = 0, hjust = 1),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 7, hjust = -0.7)) #+


matrix <- ggplot() +
  geom_density_ridges_gradient(data = Spie_posterior,
                               aes(x = S_PIE + S_PIE_global,
                                   y = Matrix.category,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(S_PIE_global, probs = 0.025),
                        upper = quantile(S_PIE_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_vline(data = frag_global,
             aes(xintercept = median(S_PIE_global)),
             size = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5) +
  geom_point(data = Spie_posterior,
             aes(x = S_PIE + S_PIE_global, 
                 y = Matrix.category),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 1.25, shape = 18) +
  geom_text(data = Spie_posterior %>%
              group_by(Matrix.category) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(Matrix.category, n_study, .keep_all = T),
            aes(x=-0.175, y=Matrix.category, 
                label=paste('n[study] == ', n_study)),
            size=1.5,
            nudge_y = 0.15, parse = T) +
  theme_bw() +
  labs(y = 'Matrix filter',
       x = '',#,#expression(paste('Study-level slope')),
       subtitle = 'Matrix filter',
       tag = 'h'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                               '#969696', '#cccccc')) +
  coord_fixed(0.5/3) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(vjust = 0, hjust = 1),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 7, hjust = -0.42)) #+

bottom2 = cowplot::plot_grid(taxa,
                             continent,
                             time,
                             matrix,
                             nrow = 2) +
  cowplot::draw_label(expression(paste('Standardised evenness ~ fragment size slope estimate')), 
                      y = 0.015, size = 7)

cowplot::plot_grid(top, bottom1, bottom2, rel_heights = c(0.025,1, 1), nrow = 3) #+

# plot for 2 column width
# set local directory to save
ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision3/figures/Ex_Dat_Fig5.pdf',
       width = 120,
       height = 210,
       units = 'mm')
