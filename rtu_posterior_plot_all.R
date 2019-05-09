# COLOUR VERSION

# beta-diversity this time

# get study-level posterior samples from model with only study-level slope variation
# code to plot each of the groups we are interested in together on one figure: 
# continent, taxa, time, biome, matrix, fragment and matrix sphere (climate?)

library(tidyverse)
library(brms)
library(ggridges)

source('~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/beta_frag_posterior_wrangle.R')

# dummy plot for creating separate legend
three_grey_legend <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = Rtu_posterior,
                               aes(x = Rtu + unique(Rtu_global_slope),
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
                    labels = c('< 5%', '5-45%',  '50%')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        # legend.justification = c(1, 1),
        legend.background = element_blank()) #+
source('~/Dropbox/1current/R_random/functions/gg_legend.R')
legend <- gg_legend(three_grey_legend)

Rtu_study_posterior_time <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = Rtu_posterior,
                               aes(x = Rtu + unique(Rtu_global_slope),
                                   y = time.since.fragmentation,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  # stat_density_ridges(data = Rtu_posterior,
  #                     geom = 'density_ridges_gradient',
  #                     aes(x = Rtu + unique(Rtu_global_slope), 
  #                     y = time.since.fragmentation,
  #                     fill = 0.5 - abs(0.5 -..ecdf..)),
  #                     calc_ecdf = T,
  #                     scale = 0.9, alpha = 0.5) +
  geom_rect(data = Rtu_posterior %>% distinct(Rtu_lower_slope, Rtu_upper_slope),
            aes(xmin = Rtu_lower_slope, xmax = Rtu_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.6) +
  geom_point(data = Rtu_posterior,
             aes(x = Rtu + unique(Rtu_global_slope), 
                 y = time.since.fragmentation),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = Rtu_posterior,
             aes(xintercept = Rtu_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Rtu_posterior %>%
              group_by(time.since.fragmentation) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(time.since.fragmentation, n_study, .keep_all = T),
            aes(x=Inf, y=time.since.fragmentation, 
                label=paste('n[study] == ', n_study)),
            size=3.5,
            nudge_y = 0.15, parse = T, hjust = 1.25) +
  theme_bw() +
  labs(y = 'Time since fragmentation',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'd'
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


Rtu_study_posterior_matrix <- ggplot() +
  geom_density_ridges_gradient(data = Rtu_posterior,
                               aes(x = Rtu + unique(Rtu_global_slope),
                                   y = Matrix.category,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = Rtu_posterior %>% distinct(Rtu_lower_slope, Rtu_upper_slope),
            aes(xmin = Rtu_lower_slope, xmax = Rtu_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.6) +
  geom_point(data = Rtu_posterior,
             aes(x = Rtu + unique(Rtu_global_slope), 
                 y = Matrix.category),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = Rtu_posterior,
             aes(xintercept = Rtu_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Rtu_posterior %>%
              group_by(Matrix.category) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(Matrix.category, n_study, .keep_all = T),
            aes(x=Inf, y=Matrix.category, 
                label=paste('n[study] == ', n_study)),
            size=3.5,
            nudge_y = 0.15, parse = T, hjust = 1.25) +
  theme_bw() +
  labs(y = 'Matrix filter',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'e'
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


Rtu_study_posterior_biome <- ggplot() +
  geom_density_ridges_gradient(data = Rtu_posterior,
                               aes(x = Rtu + unique(Rtu_global_slope),
                                   y = biome,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = Rtu_posterior %>% distinct(Rtu_lower_slope, Rtu_upper_slope),
            aes(xmin = Rtu_lower_slope, xmax = Rtu_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.6) +
  geom_point(data = Rtu_posterior,
             aes(x = Rtu + unique(Rtu_global_slope), 
                 y = biome),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = Rtu_posterior,
             aes(xintercept = Rtu_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Rtu_posterior %>%
              group_by(biome) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(biome, n_study, .keep_all = T),
            aes(x=Inf, y=biome, 
                label=paste('n[study] == ', n_study)),
            size=3.5,
            nudge_y = 0.15, parse = T, hjust = 1.25) +
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
        legend.background = element_blank()) #+

Rtu_study_posterior_taxa <- ggplot() +
  geom_density_ridges_gradient(data = Rtu_posterior,
                               aes(x = Rtu + unique(Rtu_global_slope),
                                   y = taxa,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = Rtu_posterior %>% distinct(Rtu_lower_slope, Rtu_upper_slope),
            aes(xmin = Rtu_lower_slope, xmax = Rtu_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.6) +
  geom_point(data = Rtu_posterior,
             aes(x = Rtu + unique(Rtu_global_slope), 
                 y = taxa),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = Rtu_posterior,
             aes(xintercept = Rtu_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Rtu_posterior %>%
              group_by(taxa) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(taxa, n_study, .keep_all = T),
            aes(x=Inf, y=taxa, 
                label=paste('n[study] == ', n_study)),
            size=3.5,
            nudge_y = 0.15, parse = T, hjust = 1.25) +
  theme_bw() +
  labs(y = 'Taxa',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'),
       tag = 'd'
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

Rtu_study_posterior_continent <- ggplot() +
  geom_density_ridges_gradient(data = Rtu_posterior,
                               aes(x = Rtu + unique(Rtu_global_slope),
                                   y = continent8,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = Rtu_posterior %>% distinct(Rtu_lower_slope, Rtu_upper_slope),
            aes(xmin = Rtu_lower_slope, xmax = Rtu_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.6) +
  geom_point(data = Rtu_posterior,
             aes(x = Rtu + unique(Rtu_global_slope), 
                 y = continent8),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = Rtu_posterior,
             aes(xintercept = Rtu_global_slope)) +
  geom_text(data = Rtu_posterior %>%
              group_by(continent8) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(continent8, n_study, .keep_all = T),
            aes(x=Inf, y=continent8, 
                label=paste('n[study] == ', n_study)),
            size=3.5,
            nudge_y = 0.15, parse = T, hjust = 1.25) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  labs(y = 'Continent',
       x = '',#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       tag = 'a'
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

Rtu_study_posterior_climate <- ggplot() +
  geom_density_ridges_gradient(data = Rtu_posterior,
                               aes(x = Rtu + unique(Rtu_global_slope),
                                   y = climate,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = Rtu_posterior %>% distinct(Rtu_lower_slope, Rtu_upper_slope),
            aes(xmin = Rtu_lower_slope, xmax = Rtu_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.6) +
  geom_point(data = Rtu_posterior,
             aes(x = Rtu + unique(Rtu_global_slope), 
                 y = climate),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = Rtu_posterior,
             aes(xintercept = Rtu_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Rtu_posterior %>%
              group_by(climate) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(climate, n_study, .keep_all = T),
            aes(x=Inf, y=climate, 
                label=paste('n[study] == ', n_study)),
            size=3.5,
            nudge_y = 0.15, parse = T, hjust = 1.25) +
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
        legend.background = element_blank()) #+

Rtu_study_posterior_sphere.frag <- ggplot() +
  geom_density_ridges_gradient(data = Rtu_posterior,
                               aes(x = Rtu + unique(Rtu_global_slope),
                                   y = frag_matrix,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.5,
                               linetype = 0) +
  geom_rect(data = Rtu_posterior %>% distinct(Rtu_lower_slope, Rtu_upper_slope),
            aes(xmin = Rtu_lower_slope, xmax = Rtu_upper_slope), ymin = -Inf, ymax = Inf,
            alpha = 0.6) +
  geom_point(data = Rtu_posterior,
             aes(x = Rtu + unique(Rtu_global_slope), 
                 y = frag_matrix),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 18) +
  geom_vline(data = Rtu_posterior,
             aes(xintercept = Rtu_global_slope)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = Rtu_posterior %>%
              group_by(frag_matrix) %>% 
              mutate(n_study = n_distinct(dataset_label)) %>% 
              ungroup() %>% 
              distinct(climate, n_study, .keep_all = T),
            aes(x=Inf, y=frag_matrix, 
                label=paste('n[study] == ', n_study)),
            size=3.5,
            nudge_y = 0.15, parse = T, hjust = 1.25) +
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
        legend.background = element_blank()) #+


top = cowplot::plot_grid(NULL, legend, NULL,
                         nrow = 1, rel_widths = c(0.2,0.6,0.2))
bottom = cowplot::plot_grid(Rtu_study_posterior_continent,
                            Rtu_study_posterior_biome,
                            Rtu_study_posterior_sphere.frag,
                            Rtu_study_posterior_taxa,
                            Rtu_study_posterior_time, 
                            Rtu_study_posterior_matrix, 
                            nrow = 2)

cowplot::plot_grid(top, bottom, rel_heights = c(0.05,1), nrow = 2) +
  cowplot::draw_label('Study-level balanced abundance change', y = 0.02)

ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/bal_abund_posteriors.png',
       width = 290,
       height = 220,
       units = 'mm')
