# get study-level posterior samples from model with only study-level slope variation
# code to plot more detail for each group:
# time
library(tidyverse)
library(brms)
library(ggridges)

source('~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/alpha_frag_posterior_wrangle.R')

sstd_study_posterior_tt <- ggplot() +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(ymin = Sstd_lower_slope, ymax = Sstd_upper_slope), xmin = -Inf, xmax = Inf,
            alpha = 0.75) +
  geom_point(data = Sstd_posterior %>% 
               group_by(time.since.fragmentation, taxa) %>% 
               summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                         lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                         lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                         upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                         upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
             aes(y = median_slope, x = time.since.fragmentation,
                 group = interaction(time.since.fragmentation, taxa),
                 colour = taxa),
             position = position_dodge(width = 0.35),
             size = 3) +
  # 75% CI
  geom_linerange(data = Sstd_posterior %>% 
                   group_by(time.since.fragmentation, taxa) %>% 
                   summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                             lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                             lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                             upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                             upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
                 aes(ymin = lower25, ymax = upper75,
                     x = time.since.fragmentation,
                     group = interaction(time.since.fragmentation, taxa),
                     colour = taxa),
                 position = position_dodge(width = 0.35),
                 size = 1.5) +
  # 95% CI
  geom_linerange(data = Sstd_posterior %>% 
                   group_by(time.since.fragmentation, taxa) %>% 
                   summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                             lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                             lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                             upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                             upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
                 aes(ymin = lower2.5, ymax = upper97.5,
                     x = time.since.fragmentation,
                     group = interaction(time.since.fragmentation, taxa),
                     colour = taxa),
                 position = position_dodge(width = 0.35),
                 size = 0.75) +
  geom_hline(data = Sstd_posterior,
             aes(yintercept = Sstd_global_slope)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text(data = Sstd_posterior %>%
              group_by(time.since.fragmentation, taxa) %>%
              summarise(n_study = n_distinct(dataset_label)), 
            aes(y = -0.125, x=time.since.fragmentation,
                label=paste('n == ', n_study),
                group = interaction(time.since.fragmentation, taxa),
                colour = taxa),
            size = 2.5,
            position = position_dodge(width = 0.35),
            parse = T) +
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


sstd_study_posterior_tm <- ggplot() +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(ymin = Sstd_lower_slope, ymax = Sstd_upper_slope), xmin = -Inf, xmax = Inf,
            alpha = 0.75) +
  geom_point(data = Sstd_posterior %>% 
               group_by(time.since.fragmentation, Matrix.category) %>% 
               summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                         lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                         lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                         upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                         upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
             aes(y = median_slope, x = time.since.fragmentation,
                 group = interaction(time.since.fragmentation, Matrix.category),
                 colour = Matrix.category),
             position = position_dodge(width = 0.25),
             size = 3) +
  # 75% CI
  geom_linerange(data = Sstd_posterior %>% 
                   group_by(time.since.fragmentation, Matrix.category) %>% 
                   summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                             lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                             lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                             upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                             upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
                 aes(ymin = lower25, ymax = upper75,
                     x = time.since.fragmentation,
                     group = interaction(time.since.fragmentation, Matrix.category),
                     colour = Matrix.category),
                 position = position_dodge(width = 0.25),
                 size = 1.5) +
  # 95% CI
  geom_linerange(data = Sstd_posterior %>% 
                   group_by(time.since.fragmentation, Matrix.category) %>% 
                   summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                             lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                             lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                             upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                             upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
                 aes(ymin = lower2.5, ymax = upper97.5,
                     x = time.since.fragmentation,
                     group = interaction(time.since.fragmentation, Matrix.category),
                     colour = Matrix.category),
                 position = position_dodge(width = 0.25),
                 size = 0.75) +
  geom_hline(data = Sstd_posterior,
             aes(yintercept = Sstd_global_slope)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text(data = Sstd_posterior %>%
              group_by(time.since.fragmentation, Matrix.category) %>%
              summarise(n_study = n_distinct(dataset_label)), 
            aes(y = -0.125, x=time.since.fragmentation,
                label=paste('n == ', n_study),
                group = interaction(time.since.fragmentation, Matrix.category),
                colour = Matrix.category),
            size = 2.5,
            position = position_dodge(width = 0.25),
            parse = T) +
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


sstd_study_posterior_tb <- ggplot() +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(ymin = Sstd_lower_slope, ymax = Sstd_upper_slope), xmin = -Inf, xmax = Inf,
            alpha = 0.75) +
  geom_point(data = Sstd_posterior %>% 
               group_by(time.since.fragmentation, biome) %>% 
               summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                         lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                         lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                         upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                         upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
             aes(y = median_slope, x = time.since.fragmentation,
                 group = interaction(time.since.fragmentation, biome),
                 colour = biome),
             position = position_dodge(width = 0.25),
             size = 3) +
  # 75% CI
  geom_linerange(data = Sstd_posterior %>% 
                   group_by(time.since.fragmentation, biome) %>% 
                   summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                             lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                             lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                             upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                             upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
                 aes(ymin = lower25, ymax = upper75,
                     x = time.since.fragmentation,
                     group = interaction(time.since.fragmentation, biome),
                     colour = biome),
                 position = position_dodge(width = 0.25),
                 size = 1.5) +
  # 95% CI
  geom_linerange(data = Sstd_posterior %>% 
                   group_by(time.since.fragmentation, biome) %>% 
                   summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                             lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                             lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                             upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                             upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
                 aes(ymin = lower2.5, ymax = upper97.5,
                     x = time.since.fragmentation,
                     group = interaction(time.since.fragmentation, biome),
                     colour = biome),
                 position = position_dodge(width = 0.25),
                 size = 0.75) +
  geom_hline(data = Sstd_posterior,
             aes(yintercept = Sstd_global_slope)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text(data = Sstd_posterior %>%
              group_by(time.since.fragmentation, biome) %>%
              summarise(n_study = n_distinct(dataset_label)), 
            aes(y = -0.125, x=time.since.fragmentation,
                label=paste('n == ', n_study),
                group = interaction(time.since.fragmentation, biome),
                colour = biome),
            size = 2.5,
            position = position_dodge(width = 0.25),
            parse = T) +
  theme_bw() +
  labs(x = '',
       y = ''#,#expression(paste('Study-level slope')),
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
  ) +
  scale_x_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  scale_fill_brewer(name = 'Biome', type = 'qual', palette = 'Set2') +
  scale_colour_brewer(name = 'Biome', type = 'qual', palette = 'Set2') +
  coord_flip() +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.key.size = unit(0.1, 'mm')) +
  guides(fill = guide_legend(reverse = T, nrow = 2),
         colour = guide_legend(reverse = T, nrow = 2))


sstd_study_posterior_tc <- ggplot() +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(ymin = Sstd_lower_slope, ymax = Sstd_upper_slope), xmin = -Inf, xmax = Inf,
            alpha = 0.75) +
  geom_point(data = Sstd_posterior %>% 
               group_by(time.since.fragmentation, continent8) %>% 
               summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                         lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                         lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                         upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                         upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
             aes(y = median_slope, x = time.since.fragmentation,
                 group = interaction(time.since.fragmentation, continent8),
                 colour = continent8),
             position = position_dodge(width = 0.5),
             size = 3) +
  # 75% CI
  geom_linerange(data = Sstd_posterior %>% 
                   group_by(time.since.fragmentation, continent8) %>% 
                   summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                             lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                             lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                             upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                             upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
                 aes(ymin = lower25, ymax = upper75,
                     x = time.since.fragmentation,
                     group = interaction(time.since.fragmentation, continent8),
                     colour = continent8),
                 position = position_dodge(width = 0.5),
                 size = 1.5) +
  # 95% CI
  geom_linerange(data = Sstd_posterior %>% 
                   group_by(time.since.fragmentation, continent8) %>% 
                   summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                             lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                             lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                             upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                             upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
                 aes(ymin = lower2.5, ymax = upper97.5,
                     x = time.since.fragmentation,
                     group = interaction(time.since.fragmentation, continent8),
                     colour = continent8),
                 position = position_dodge(width = 0.5),
                 size = 0.75) +
  geom_hline(data = Sstd_posterior,
             aes(yintercept = Sstd_global_slope)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text(data = Sstd_posterior %>%
              group_by(time.since.fragmentation, continent8) %>%
              summarise(n_study = n_distinct(dataset_label)), 
            aes(y = -0.125, x=time.since.fragmentation,
                label=paste('n == ', n_study),
                group = interaction(time.since.fragmentation, continent8),
                colour = continent8),
            size = 2.5,
            position = position_dodge(width = 0.5),
            parse = T) +
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


sstd_study_posterior_tl <- ggplot() +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(ymin = Sstd_lower_slope, ymax = Sstd_upper_slope), xmin = -Inf, xmax = Inf,
            alpha = 0.75) +
  geom_point(data = Sstd_posterior %>% 
               group_by(time.since.fragmentation, climate) %>% 
               summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                         lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                         lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                         upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                         upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
             aes(y = median_slope, x = time.since.fragmentation,
                 group = interaction(time.since.fragmentation, climate),
                 colour = climate),
             position = position_dodge(width = 0.25),
             size = 3) +
  # 75% CI
  geom_linerange(data = Sstd_posterior %>% 
                   group_by(time.since.fragmentation, climate) %>% 
                   summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                             lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                             lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                             upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                             upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
                 aes(ymin = lower25, ymax = upper75,
                     x = time.since.fragmentation,
                     group = interaction(time.since.fragmentation, climate),
                     colour = climate),
                 position = position_dodge(width = 0.25),
                 size = 1.5) +
  # 95% CI
  geom_linerange(data = Sstd_posterior %>% 
                   group_by(time.since.fragmentation, climate) %>% 
                   summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                             lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                             lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                             upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                             upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
                 aes(ymin = lower2.5, ymax = upper97.5,
                     x = time.since.fragmentation,
                     group = interaction(time.since.fragmentation, climate),
                     colour = climate),
                 position = position_dodge(width = 0.25),
                 size = 0.75) +
  geom_hline(data = Sstd_posterior,
             aes(yintercept = Sstd_global_slope)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text(data = Sstd_posterior %>%
              group_by(time.since.fragmentation, climate) %>%
              summarise(n_study = n_distinct(dataset_label)), 
            aes(y = -0.125, x=time.since.fragmentation,
                label=paste('n == ', n_study),
                group = interaction(time.since.fragmentation, climate),
                colour = climate),
            size = 2.5,
            position = position_dodge(width = 0.25),
            parse = T) +
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

sstd_study_posterior_ts <- ggplot() +
  geom_rect(data = Sstd_posterior %>% distinct(Sstd_lower_slope, Sstd_upper_slope),
            aes(ymin = Sstd_lower_slope, ymax = Sstd_upper_slope), xmin = -Inf, xmax = Inf,
            alpha = 0.75) +
  geom_point(data = Sstd_posterior %>% 
               group_by(time.since.fragmentation, frag_matrix) %>% 
               summarise(median_slope = median(S_std + unique(Sstd_global_slope))),
             aes(y = median_slope, x = time.since.fragmentation,
                 group = interaction(time.since.fragmentation, frag_matrix),
                 colour = frag_matrix),
             position = position_dodge(width = 0.25),
             size = 3) +
  # 75% CI
  geom_linerange(data = Sstd_posterior %>% 
                   group_by(time.since.fragmentation, frag_matrix) %>% 
                   summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                             lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                             lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                             upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                             upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
                 aes(ymin = lower25, ymax = upper75,
                     x = time.since.fragmentation,
                     group = interaction(time.since.fragmentation, frag_matrix),
                     colour = frag_matrix),
                 position = position_dodge(width = 0.25),
                 size = 1.5) +
  # 95% CI
  geom_linerange(data = Sstd_posterior %>% 
                   group_by(time.since.fragmentation, frag_matrix) %>% 
                   summarise(median_slope = median(S_std + unique(Sstd_global_slope)),
                             lower2.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.025),
                             lower25 = quantile(S_std + unique(Sstd_global_slope), probs = 0.25),
                             upper75 = quantile(S_std + unique(Sstd_global_slope), probs = 0.75),
                             upper97.5 = quantile(S_std + unique(Sstd_global_slope), probs = 0.975)),
                 aes(ymin = lower2.5, ymax = upper97.5,
                     x = time.since.fragmentation,
                     group = interaction(time.since.fragmentation, frag_matrix),
                     colour = frag_matrix),
                 position = position_dodge(width = 0.25),
                 size = 0.75) +
  geom_hline(data = Sstd_posterior,
             aes(yintercept = Sstd_global_slope)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text(data = Sstd_posterior %>%
              group_by(time.since.fragmentation, frag_matrix) %>%
              summarise(n_study = n_distinct(dataset_label)), 
            aes(y = -0.125, x=time.since.fragmentation,
                label=paste('n == ', n_study),
                group = interaction(time.since.fragmentation, frag_matrix),
                colour = frag_matrix),
            size = 2.5,
            position = position_dodge(width = 0.25),
            parse = T) +
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


bottom <- cowplot::plot_grid(sstd_study_posterior_tl,
                             sstd_study_posterior_tc,
                             sstd_study_posterior_tb,
                             sstd_study_posterior_ts,
                             sstd_study_posterior_tt,
                             sstd_study_posterior_tm
                             )
bottom + 
  cowplot::draw_label('Study-level slope', y = 0.01) +
  cowplot::draw_label('Time since fragmentation', x = 0.01, angle = 90)

ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/FigSx_time.png',
       width = 350,
       height = 220,
       units = 'mm')
