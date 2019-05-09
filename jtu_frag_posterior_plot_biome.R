# get study-level posterior samples from model with only study-level slope variation
# code to plot more detail for each group:
# biome
library(tidyverse)
library(brms)
library(ggridges)

source('~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/beta_frag_posterior_wrangle.R')

Jtu_study_posterior_tt <- ggplot() +
  geom_rect(data = Jtu_posterior %>% distinct(Jtu_lower_slope, Jtu_upper_slope),
            aes(ymin = Jtu_lower_slope, ymax = Jtu_upper_slope), xmin = -Inf, xmax = Inf,
            alpha = 0.75) +
  geom_point(data = Jtu_posterior %>% 
               group_by(biome, taxa) %>% 
               summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                         lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                         lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                         upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                         upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
             aes(y = median_slope, x = biome,
                 group = interaction(biome, taxa),
                 colour = taxa),
             position = position_dodge(width = 0.85),
             size = 3) +
  # 75% CI
  geom_linerange(data = Jtu_posterior %>% 
                   group_by(biome, taxa) %>% 
                   summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                             lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                             lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                             upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                             upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
                 aes(ymin = lower25, ymax = upper75,
                     x = biome,
                     group = interaction(biome, taxa),
                     colour = taxa),
                 position = position_dodge(width = 0.85),
                 size = 1.5) +
  # 95% CI
  geom_linerange(data = Jtu_posterior %>% 
                   group_by(biome, taxa) %>% 
                   summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                             lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                             lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                             upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                             upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
                 aes(ymin = lower2.5, ymax = upper97.5,
                     x = biome,
                     group = interaction(biome, taxa),
                     colour = taxa),
                 position = position_dodge(width = 0.85),
                 size = 0.75) +
  geom_hline(data = Jtu_posterior,
             aes(yintercept = Jtu_global_slope)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text(data = Jtu_posterior %>%
              group_by(biome, taxa) %>%
              summarise(n_study = n_distinct(dataset_label)), 
            aes(y = Inf, x=biome,
                label=paste('n == ', n_study),
                group = interaction(biome, taxa),
                colour = taxa),
            size = 2.5,
            position = position_dodge(width = 0.85),
            hjust = 1.25, parse = T) +
  theme_bw() +
  labs(x = '',
       y = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_x_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  scale_fill_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  coord_flip() +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.key.size = unit(0.1, 'mm')) +
  guides(fill = guide_legend(reverse = T, nrow = 2),
         colour = guide_legend(reverse = T, nrow = 2))


Jtu_study_posterior_tm <- ggplot() +
  geom_rect(data = Jtu_posterior %>% distinct(Jtu_lower_slope, Jtu_upper_slope),
            aes(ymin = Jtu_lower_slope, ymax = Jtu_upper_slope), xmin = -Inf, xmax = Inf,
            alpha = 0.75) +
  geom_point(data = Jtu_posterior %>% 
               group_by(biome, Matrix.category) %>% 
               summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                         lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                         lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                         upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                         upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
             aes(y = median_slope, x = biome,
                 group = interaction(biome, Matrix.category),
                 colour = Matrix.category),
             position = position_dodge(width = 0.5),
             size = 3) +
  # 75% CI
  geom_linerange(data = Jtu_posterior %>% 
                   group_by(biome, Matrix.category) %>% 
                   summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                             lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                             lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                             upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                             upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
                 aes(ymin = lower25, ymax = upper75,
                     x = biome,
                     group = interaction(biome, Matrix.category),
                     colour = Matrix.category),
                 position = position_dodge(width = 0.5),
                 size = 1.5) +
  # 95% CI
  geom_linerange(data = Jtu_posterior %>% 
                   group_by(biome, Matrix.category) %>% 
                   summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                             lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                             lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                             upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                             upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
                 aes(ymin = lower2.5, ymax = upper97.5,
                     x = biome,
                     group = interaction(biome, Matrix.category),
                     colour = Matrix.category),
                 position = position_dodge(width = 0.5),
                 size = 0.75) +
  geom_hline(data = Jtu_posterior,
             aes(yintercept = Jtu_global_slope)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text(data = Jtu_posterior %>%
              group_by(biome, Matrix.category) %>%
              summarise(n_study = n_distinct(dataset_label)), 
            aes(y = Inf, x=biome,
                label=paste('n == ', n_study),
                group = interaction(biome, Matrix.category),
                colour = Matrix.category),
            size = 2.5,
            position = position_dodge(width = 0.5),
            hjust = 1.25, parse = T) +
  theme_bw() +
  labs(x = '',
       y = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_x_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  scale_fill_brewer(name = 'Matrix filter', type = 'qual', palette = 'Accent') +
  scale_colour_brewer(name = 'Matrix filter', type = 'qual', palette = 'Accent') +
  coord_flip() +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.key.size = unit(0.1, 'mm')) +
  guides(fill = guide_legend(reverse = T, nrow = 2),
         colour = guide_legend(reverse = T, nrow = 2))


Jtu_study_posterior_tb <- ggplot() +
  geom_rect(data = Jtu_posterior %>% distinct(Jtu_lower_slope, Jtu_upper_slope),
            aes(ymin = Jtu_lower_slope, ymax = Jtu_upper_slope), xmin = -Inf, xmax = Inf,
            alpha = 0.75) +
  geom_point(data = Jtu_posterior %>% 
               group_by(biome, continent8) %>% 
               summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                         lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                         lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                         upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                         upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
             aes(y = median_slope, x = biome,
                 group = interaction(biome, continent8),
                 colour = continent8),
             position = position_dodge(width = 0.8),
             size = 3) +
  # 75% CI
  geom_linerange(data = Jtu_posterior %>% 
                   group_by(biome, continent8) %>% 
                   summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                             lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                             lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                             upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                             upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
                 aes(ymin = lower25, ymax = upper75,
                     x = biome,
                     group = interaction(biome, continent8),
                     colour = continent8),
                 position = position_dodge(width = 0.8),
                 size = 1.5) +
  # 95% CI
  geom_linerange(data = Jtu_posterior %>% 
                   group_by(biome, continent8) %>% 
                   summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                             lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                             lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                             upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                             upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
                 aes(ymin = lower2.5, ymax = upper97.5,
                     x = biome,
                     group = interaction(biome, continent8),
                     colour = continent8),
                 position = position_dodge(width = 0.8),
                 size = 0.75) +
  geom_hline(data = Jtu_posterior,
             aes(yintercept = Jtu_global_slope)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text(data = Jtu_posterior %>%
              group_by(biome, continent8) %>%
              summarise(n_study = n_distinct(dataset_label)), 
            aes(y = Inf, x=biome,
                label=paste('n == ', n_study),
                group = interaction(biome, continent8),
                colour = continent8),
            size = 2.5,
            position = position_dodge(width = 0.8),
            hjust = 1.25, parse = T) +
  theme_bw() +
  labs(x = '',
       y = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_x_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  scale_fill_brewer(name = 'Region', type = 'qual', palette = 'Paired',
                    breaks = c('Africa', 'Asia', 'Central America',
                               'Europe', 'North America', 'Oceania',
                               'South America'),
                    labels = c('Africa', 'Asia',
                               'C America', 'Europe', 'N America',
                               'Oceania', 'S America')) +
  scale_colour_brewer(name = 'Region', type = 'qual', palette = 'Paired',
                      breaks = c('Africa', 'Asia', 'Central America',
                                 'Europe', 'North America', 'Oceania',
                                 'South America'),
                      labels = c('Africa', 'Asia',
                                 'C America', 'Europe', 'N America',
                                 'Oceania', 'S America')) +
  coord_flip() +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.key.size = unit(0.1, 'mm')) +
  guides(fill = guide_legend(reverse = T, nrow = 2),
         colour = guide_legend(reverse = T, nrow = 2))


Jtu_study_posterior_tc <- ggplot() +
  geom_rect(data = Jtu_posterior %>% distinct(Jtu_lower_slope, Jtu_upper_slope),
            aes(ymin = Jtu_lower_slope, ymax = Jtu_upper_slope), xmin = -Inf, xmax = Inf,
            alpha = 0.75) +
  geom_point(data = Jtu_posterior %>% 
               group_by(biome, time.since.fragmentation) %>% 
               summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                         lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                         lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                         upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                         upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
             aes(y = median_slope, x = biome,
                 group = interaction(biome, time.since.fragmentation),
                 colour = time.since.fragmentation),
             position = position_dodge(width = 0.5),
             size = 3) +
  # 75% CI
  geom_linerange(data = Jtu_posterior %>% 
                   group_by(biome, time.since.fragmentation) %>% 
                   summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                             lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                             lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                             upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                             upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
                 aes(ymin = lower25, ymax = upper75,
                     x = biome,
                     group = interaction(biome, time.since.fragmentation),
                     colour = time.since.fragmentation),
                 position = position_dodge(width = 0.5),
                 size = 1.5) +
  # 95% CI
  geom_linerange(data = Jtu_posterior %>% 
                   group_by(biome, time.since.fragmentation) %>% 
                   summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                             lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                             lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                             upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                             upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
                 aes(ymin = lower2.5, ymax = upper97.5,
                     x = biome,
                     group = interaction(biome, time.since.fragmentation),
                     colour = time.since.fragmentation),
                 position = position_dodge(width = 0.5),
                 size = 0.75) +
  geom_hline(data = Jtu_posterior,
             aes(yintercept = Jtu_global_slope)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text(data = Jtu_posterior %>%
              group_by(biome, time.since.fragmentation) %>%
              summarise(n_study = n_distinct(dataset_label)), 
            aes(y = Inf, x=biome,
                label=paste('n == ', n_study),
                group = interaction(biome, time.since.fragmentation),
                colour = time.since.fragmentation),
            size = 2.5,
            position = position_dodge(width = 0.5),
            hjust = 1.25, parse = T) +
  theme_bw() +
  labs(x = '',
       y = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_x_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  scale_fill_brewer(name = 'Time since\nfragmentation', type = 'qual', palette = 'Paired') +
  scale_colour_brewer(name = 'Time since\nfragmentation', type = 'qual', palette = 'Paired') +
  coord_flip() +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.key.size = unit(0.1, 'mm')) +
  guides(fill = guide_legend(reverse = T, nrow = 2),
         colour = guide_legend(reverse = T, nrow = 2))


Jtu_study_posterior_tl <- ggplot() +
  geom_rect(data = Jtu_posterior %>% distinct(Jtu_lower_slope, Jtu_upper_slope),
            aes(ymin = Jtu_lower_slope, ymax = Jtu_upper_slope), xmin = -Inf, xmax = Inf,
            alpha = 0.75) +
  geom_point(data = Jtu_posterior %>% 
               group_by(biome, climate) %>% 
               summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                         lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                         lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                         upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                         upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
             aes(y = median_slope, x = biome,
                 group = interaction(biome, climate),
                 colour = climate),
             position = position_dodge(width = 0.5),
             size = 3) +
  # 75% CI
  geom_linerange(data = Jtu_posterior %>% 
                   group_by(biome, climate) %>% 
                   summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                             lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                             lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                             upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                             upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
                 aes(ymin = lower25, ymax = upper75,
                     x = biome,
                     group = interaction(biome, climate),
                     colour = climate),
                 position = position_dodge(width = 0.5),
                 size = 1.5) +
  # 95% CI
  geom_linerange(data = Jtu_posterior %>% 
                   group_by(biome, climate) %>% 
                   summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                             lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                             lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                             upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                             upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
                 aes(ymin = lower2.5, ymax = upper97.5,
                     x = biome,
                     group = interaction(biome, climate),
                     colour = climate),
                 position = position_dodge(width = 0.5),
                 size = 0.75) +
  geom_hline(data = Jtu_posterior,
             aes(yintercept = Jtu_global_slope)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text(data = Jtu_posterior %>%
              group_by(biome, climate) %>%
              summarise(n_study = n_distinct(dataset_label)), 
            aes(y = Inf, x=biome,
                label=paste('n == ', n_study),
                group = interaction(biome, climate),
                colour = climate),
            size = 2.5,
            position = position_dodge(width = 0.5),
            hjust = 1.25, parse = T) +
  theme_bw() +
  labs(x = '',
       y = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_x_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  scale_fill_brewer(name = 'Climate', type = 'qual', palette = 'Accent') +
  scale_colour_brewer(name = 'Climate', type = 'qual', palette = 'Accent') +
  coord_flip() +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.key.size = unit(0.1, 'mm')) +
  guides(fill = guide_legend(reverse = T, nrow = 2),
         colour = guide_legend(reverse = T, nrow = 2))

Jtu_study_posterior_ts <- ggplot() +
  geom_rect(data = Jtu_posterior %>% distinct(Jtu_lower_slope, Jtu_upper_slope),
            aes(ymin = Jtu_lower_slope, ymax = Jtu_upper_slope), xmin = -Inf, xmax = Inf,
            alpha = 0.75) +
  geom_point(data = Jtu_posterior %>% 
               group_by(biome, frag_matrix) %>% 
               summarise(median_slope = median(Jtu + unique(Jtu_global_slope))),
             aes(y = median_slope, x = biome,
                 group = interaction(biome, frag_matrix),
                 colour = frag_matrix),
             position = position_dodge(width = 0.5),
             size = 3) +
  # 75% CI
  geom_linerange(data = Jtu_posterior %>% 
                   group_by(biome, frag_matrix) %>% 
                   summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                             lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                             lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                             upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                             upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
                 aes(ymin = lower25, ymax = upper75,
                     x = biome,
                     group = interaction(biome, frag_matrix),
                     colour = frag_matrix),
                 position = position_dodge(width = 0.5),
                 size = 1.5) +
  # 95% CI
  geom_linerange(data = Jtu_posterior %>% 
                   group_by(biome, frag_matrix) %>% 
                   summarise(median_slope = median(Jtu + unique(Jtu_global_slope)),
                             lower2.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.025),
                             lower25 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.25),
                             upper75 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.75),
                             upper97.5 = quantile(Jtu + unique(Jtu_global_slope), probs = 0.975)),
                 aes(ymin = lower2.5, ymax = upper97.5,
                     x = biome,
                     group = interaction(biome, frag_matrix),
                     colour = frag_matrix),
                 position = position_dodge(width = 0.5),
                 size = 0.75) +
  geom_hline(data = Jtu_posterior,
             aes(yintercept = Jtu_global_slope)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text(data = Jtu_posterior %>%
              group_by(biome, frag_matrix) %>%
              summarise(n_study = n_distinct(dataset_label)), 
            aes(y = Inf, x=biome,
                label=paste('n == ', n_study),
                group = interaction(biome, frag_matrix),
                colour = frag_matrix),
            size = 2.5,
            position = position_dodge(width = 0.5),
            hjust = 1.25, parse = T) +
  theme_bw() +
  labs(x = '',
       y = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_x_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  scale_fill_brewer(name = 'Fragment &\nmatrix sphere', type = 'qual', palette = 'Accent') +
  scale_colour_brewer(name = 'Fragment &\nmatrix sphere', type = 'qual', palette = 'Accent') +
  coord_flip() +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.key.size = unit(0.1, 'mm')) +
  guides(fill = guide_legend(reverse = T, nrow = 2),
         colour = guide_legend(reverse = T, nrow = 2))


bottom <- cowplot::plot_grid(Jtu_study_posterior_tl,
                             Jtu_study_posterior_tb,
                             Jtu_study_posterior_ts,
                             Jtu_study_posterior_tt,
                             Jtu_study_posterior_tc,                             
                             Jtu_study_posterior_tm
)

bottom + 
  cowplot::draw_label('Study-level turnover', y = 0.01) +
  cowplot::draw_label('Biome', x = 0.01, angle = 90)

ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/FigSx_biome_turnover.png',
       width = 350,
       height = 220,
       units = 'mm')
