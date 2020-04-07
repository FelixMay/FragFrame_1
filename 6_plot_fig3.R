# need to execute 0_init_dirs_load_packages.R first

# COLOUR VERSION

# get study-level posterior samples from model with only study-level slope variation
# code to plot each of the groups we are interested in together on one figure: 
# taxa, continent, time, matrix

library(ggridges)

source(paste0(path2wd, '5aa_alpha_frag_posterior_wrangle.R'))

## dummy plot for creating separate legend
three_grey_legend <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = Sstd2_posterior,
                               aes(x = S_std + Sstd2_global,
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
        legend.text = element_text(size = 6, face = 'plain'),
        legend.title = element_text(size = 7, face = 'plain'),
        legend.margin = margin(),
        legend.box.spacing = unit(c(0,0,0,0), units = 'mm')) #+

source('~/Dropbox/1current/R_random/functions/gg_legend.R')
legend <- gg_legend(three_grey_legend)

Sstd2_study_posterior_time <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = Sstd2_posterior,
                               aes(x = S_std + Sstd2_global,
                                   y = time.since.fragmentation,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  # stat_density_ridges(data = Sstd2_posterior,
  #                     geom = 'density_ridges_gradient',
  #                     aes(x = S_std + unique(Sstd2_global_slope), 
  #                     y = time.since.fragmentation,
  #                     fill = 0.5 - abs(0.5 -..ecdf..)),
  #                     calc_ecdf = T,
  #                     scale = 0.9, alpha = 0.5) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Sstd2_global, probs = 0.025),
                        upper = quantile(Sstd2_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_point(data = Sstd2_posterior,
             aes(x = S_std + Sstd2_global, 
                 y = time.since.fragmentation),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Sstd2_global))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Sstd2_posterior %>%
              group_by(time.since.fragmentation) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(time.since.fragmentation, n_study, .keep_all = T),
            aes(x=0.3, y=time.since.fragmentation, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Time since fragmentation',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'c'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#bdd7e7', '#6baed6', '#3182bd',
                               '#6baed6', '#bdd7e7')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        text = element_text(size = 7),
        plot.tag = element_text(size = 8, face = 'bold')) #+


Sstd2_study_posterior_matrix <- ggplot() +
  geom_density_ridges_gradient(data = Sstd2_posterior,
                               aes(x = S_std + Sstd2_global,
                                   y = Matrix.category,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Sstd2_global, probs = 0.025),
                        upper = quantile(Sstd2_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_point(data = Sstd2_posterior,
             aes(x = S_std + Sstd2_global, 
                 y = Matrix.category),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Sstd2_global))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Sstd2_posterior %>%
              group_by(Matrix.category) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(Matrix.category, n_study, .keep_all = T),
            aes(x=0.3, y=Matrix.category, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Matrix filter',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'd'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#bae4b3', '#74c476', '#31a354',
                               '#74c476', '#bae4b3')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        text = element_text(size = 7),
        plot.tag = element_text(size = 8, face = 'bold')) #+


Sstd2_study_posterior_biome <- ggplot() +
  geom_density_ridges_gradient(data = Sstd2_posterior,
                               aes(x = S_std + Sstd2_global,
                                   y = biome,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Sstd2_global, probs = 0.025),
                        upper = quantile(Sstd2_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_point(data = Sstd2_posterior,
             aes(x = S_std + Sstd2_global, 
                 y = biome),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Sstd2_global))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Sstd2_posterior %>%
              group_by(biome) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(biome, n_study, .keep_all = T),
            aes(x=0.3, y=biome, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Biome',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'b'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#fdbe85', '#fd8d3c', '#e6550d',
                               '#fd8d3c', '#fdbe85')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        text = element_text(size = 7),
        plot.tag = element_text(size = 8, face = 'bold')) #+

Sstd2_study_posterior_taxa <- ggplot() +
  geom_density_ridges_gradient(data = Sstd2_posterior,
                               aes(x = S_std + Sstd2_global,
                                   y = taxa,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Sstd2_global, probs = 0.025),
                        upper = quantile(Sstd2_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_point(data = Sstd2_posterior,
             aes(x = S_std + Sstd2_global, 
                 y = taxa),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Sstd2_global))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Sstd2_posterior %>%
              group_by(taxa) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(taxa, n_study, .keep_all = T),
            aes(x=0.3, y=taxa, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Taxon group',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'),
       tag = 'a'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#fcae91', '#fb6a4a', '#de2d26',
                               '#fb6a4a', '#fcae91')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        text = element_text(size = 7),
        plot.tag = element_text(size = 8, face = 'bold')) #+

Sstd2_study_posterior_continent <- ggplot() +
  geom_density_ridges_gradient(data = Sstd2_posterior,
                               aes(x = S_std + Sstd2_global,
                                   y = continent8,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Sstd2_global, probs = 0.025),
                        upper = quantile(Sstd2_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_point(data = Sstd2_posterior,
             aes(x = S_std + Sstd2_global, 
                 y = continent8),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Sstd2_global))) +
  geom_text(data = Sstd2_posterior %>%
              group_by(continent8) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(continent8, n_study, .keep_all = T),
            aes(x=0.3, y=continent8, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Region',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'b'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cbc9e2', '#9e9ac8', '#756bb1',
                               '#9e9ac8', '#cbc9e2')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none',
        # legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        text = element_text(size = 7),
        plot.tag = element_text(size = 8, face = 'bold')) #+

Sstd2_study_posterior_climate <- ggplot() +
  geom_density_ridges_gradient(data = Sstd2_posterior,
                               aes(x = S_std + Sstd2_global,
                                   y = climate,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Sstd2_global, probs = 0.025),
                        upper = quantile(Sstd2_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_point(data = Sstd2_posterior,
             aes(x = S_std + Sstd2_global, 
                 y = climate),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Sstd2_global))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Sstd2_posterior %>%
              group_by(climate) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(climate, n_study, .keep_all = T),
            aes(x=0.3, y=climate, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Climate',
       x = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#636363',
                               '#969696', '#cccccc')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none',
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        text = element_text(size = 7),
        plot.tag = element_text(size = 8, face = 'bold')) #+

Sstd2_study_posterior_sphere.frag <- ggplot() +
  geom_density_ridges_gradient(data = Sstd2_posterior,
                               aes(x = S_std + Sstd2_global,
                                   y = frag_matrix,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = frag_global %>% 
              summarise(lower = quantile(Sstd2_global, probs = 0.025),
                        upper = quantile(Sstd2_global, probs = 0.975)),
            aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
            alpha = 0.6) +
  geom_point(data = Sstd2_posterior,
             aes(x = S_std + Sstd2_global, 
                 y = frag_matrix),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Sstd2_global))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Sstd2_posterior %>%
              group_by(frag_matrix) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(climate, n_study, .keep_all = T),
            aes(x=0.3, y=frag_matrix, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Fragment and matrix sphere',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'c'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#ecd1e8', '#dba7cd', '#cc7db0',
                               '#dba7cd', '#ecd1e8')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none',
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        text = element_text(size = 7),
        plot.tag = element_text(size = 8, face = 'bold')) #+


top = cowplot::plot_grid(NULL, legend, NULL,
                         nrow = 1, rel_widths = c(0.2,0.6,0.2))
bottom = cowplot::plot_grid(
                            #Sstd2_study_posterior_biome,
                            #Sstd2_study_posterior_sphere.frag,
                            Sstd2_study_posterior_taxa,
                            Sstd2_study_posterior_continent,
                            Sstd2_study_posterior_time, 
                            Sstd2_study_posterior_matrix, 
                            nrow = 2)
cowplot::plot_grid(top, bottom, rel_heights = c(0.05,1), nrow = 2) +
  cowplot::draw_label(expression(paste('Standardised species richness ~ fragment size slope estimate')), y = 0.01, size = 7)

# plot for 2 column width 
ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision3/figures/fig3_2column.pdf',
       width = 183,
       height = 170,
       units = 'mm')
# 
##repeat for Nstd and Nstd for supplement
continent <- ggplot() +
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
  geom_point(data = Nstd_posterior,
             aes(x = Nstd + Nstd_global_slope, 
                 y = continent8),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Nstd_global))) +
  geom_text(data = Nstd_posterior %>%
              group_by(continent8) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(continent8, n_study, .keep_all = T),
            aes(x=-0.2, y=continent8, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Continent',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'b'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cbc9e2', '#9e9ac8', '#756bb1',
                               '#9e9ac8', '#cbc9e2')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none',
        # legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

taxa <- ggplot() +
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
  geom_point(data = Nstd_posterior,
             aes(x = Nstd + Nstd_global_slope, 
                 y = taxa),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Nstd_global))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Nstd_posterior %>%
              group_by(taxa) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(taxa, n_study, .keep_all = T),
            aes(x=-0.2, y=taxa, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Taxon group',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'),
       tag = 'a'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#fcae91', '#fb6a4a', '#de2d26',
                               '#fb6a4a', '#fcae91')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

time <- ggplot() +
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
  geom_point(data = Nstd_posterior,
             aes(x = Nstd + Nstd_global_slope, 
                 y = time.since.fragmentation),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Nstd_global))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Nstd_posterior %>%
              group_by(time.since.fragmentation) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(time.since.fragmentation, n_study, .keep_all = T),
            aes(x=-0.2, y=time.since.fragmentation, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Time since fragmentation',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'c'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#bdd7e7', '#6baed6', '#3182bd',
                               '#6baed6', '#bdd7e7')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+


matrix <- ggplot() +
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
  geom_point(data = Nstd_posterior,
             aes(x = Nstd + Nstd_global_slope, 
                 y = Matrix.category),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(Nstd_global))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Nstd_posterior %>%
              group_by(Matrix.category) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(Matrix.category, n_study, .keep_all = T),
            aes(x=-0.2, y=Matrix.category, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Matrix filter',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'd'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#bae4b3', '#74c476', '#31a354',
                               '#74c476', '#bae4b3')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

bottom = cowplot::plot_grid(
                            #Sstd2_study_posterior_biome,
                            #Sstd2_study_posterior_sphere.frag,
                            taxa,
                            continent,
                            time,  
                            matrix,
                            nrow = 2)

cowplot::plot_grid(top, bottom, rel_heights = c(0.05,1), nrow = 2) +
  cowplot::draw_label(expression(paste('Standardised number of individuals ~ fragment size slope estimate')), y = 0.01)

# ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision1/figures/figS3_revision.png',
#        width = 240,
#        height = 220,
#        units = 'mm')

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
  geom_point(data = Spie_posterior,
             aes(x = S_PIE + S_PIE_global, 
                 y = continent8),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(S_PIE_global))) +
  geom_text(data = Spie_posterior %>%
              group_by(continent8) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(continent8, n_study, .keep_all = T),
            aes(x=-0.2, y=continent8, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Continent',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'b'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cbc9e2', '#9e9ac8', '#756bb1',
                               '#9e9ac8', '#cbc9e2')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none',
        # legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

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
  geom_point(data = Spie_posterior,
             aes(x = S_PIE + S_PIE_global, 
                 y = taxa),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(S_PIE_global))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Spie_posterior %>%
              group_by(taxa) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(taxa, n_study, .keep_all = T),
            aes(x=-0.2, y=taxa, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Taxon group',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'),
       tag = 'a'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#fcae91', '#fb6a4a', '#de2d26',
                               '#fb6a4a', '#fcae91')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

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
  geom_point(data = Spie_posterior,
             aes(x = S_PIE + S_PIE_global, 
                 y = time.since.fragmentation),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(S_PIE_global))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Spie_posterior %>%
              group_by(time.since.fragmentation) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(time.since.fragmentation, n_study, .keep_all = T),
            aes(x=-0.2, y=time.since.fragmentation, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Time since fragmentation',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'c'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#bdd7e7', '#6baed6', '#3182bd',
                               '#6baed6', '#bdd7e7')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+


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
  geom_point(data = Spie_posterior,
             aes(x = S_PIE + S_PIE_global, 
                 y = Matrix.category),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = frag_global,
             aes(xintercept = median(S_PIE_global))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Spie_posterior %>%
              group_by(Matrix.category) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(Matrix.category, n_study, .keep_all = T),
            aes(x=-0.2, y=Matrix.category, 
                label=paste('n[study] == ', n_study)),
            size=2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Matrix filter',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'd'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  # scale_fill_viridis_c(name = 'Posterior probability') +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#bae4b3', '#74c476', '#31a354',
                               '#74c476', '#bae4b3')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none', 
        legend.justification = c(1, 1),
        legend.background = element_blank()) #+

bottom = cowplot::plot_grid(
  #Sstd2_study_posterior_biome,
  #Sstd2_study_posterior_sphere.frag,
  taxa,
  continent,
  time,  
  matrix,
  nrow = 2)

cowplot::plot_grid(top, bottom, rel_heights = c(0.05,1), nrow = 2) +
  cowplot::draw_label(expression(paste('Standardised evenness ~ fragment size slope estimate')), y = 0.01)

# ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision1/figures/figS4_revision.png',
#        width = 240,
#        height = 220,
#        units = 'mm')
