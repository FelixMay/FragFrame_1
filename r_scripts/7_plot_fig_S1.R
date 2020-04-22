## code to plot simulation results for extended data

# load results of simulations
simDat <- read.csv(paste0(path2wd, 'intermediate_results/7_resultsS2000_N40000_mp1_nrep2000.csv')) %>% 
  as_tibble()


# fit linear models to each metric (for each level of aggregation)
lm_models <- simDat %>% 
  # want to fit models on a log-log scale
  mutate(logA = log(patch_area),
         logValue = log(value)) %>% 
  # throw infinite values (e.g., when value = 0)
  filter(is.finite(logValue)) %>% 
  group_by(rep, sigma, metric) %>%
  nest(logA, logValue) %>% 
  mutate(lm = purrr::map(data, ~lm(.x$logValue ~ .x$logA))) %>% 
  ungroup()

slope_coefs <- lm_models %>% 
  mutate(slope = purrr::map(lm, ~coef(.x)[2])) %>% 
  unnest(slope) %>% 
  mutate(aggregation = ifelse(sigma==1, 'Random',
                              ifelse(sigma==0.1, 'Intermediate aggregation', 'High aggregation')),
         metric2 = ifelse(metric=='S_PIE', 'Evenness',
                          ifelse(metric=='N', 'Number of individuals', 'Species richness')))

slope_coefs$aggregation <- factor(slope_coefs$aggregation,
                                  levels = c('Random', 'Intermediate aggregation', 'High aggregation'))

slope_coefs$metric2 <- factor(slope_coefs$metric2,
                              levels = c('Number of individuals', 'Species richness', 'Evenness'))

## dummy plot for creating separate legend
three_grey_legend <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = slope_coefs,
                               aes(x = slope,
                                   y = 1,
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
  scale_fill_manual(name = 'Percentiles',
                    values = c('#e5f5e0', '#a1d99b', '#31a354'),
                    labels = c('< 5%', '< 45%',  '50%')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.text = element_text(size = 6, face = 'plain'),
        legend.title = element_text(size = 7, face = 'plain'),
        legend.background = element_blank()) #+

source(paste0(path2wd, 'r_scripts/99_gg_legend.R'))
legend <- gg_legend(three_grey_legend)


sim_slopes <- ggplot() +
  facet_grid(metric2 ~ aggregation) +
  geom_density_ridges_gradient(data = slope_coefs,
                               aes(x = slope,
                                   y = 1, 
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, alpha = 0.25,
                               linetype = 0) +
  geom_point(data = slope_coefs,
             aes(x = slope, 
                 y = 1),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 3, shape = 18) +
  geom_vline(xintercept = 0, lty = 2) +
  scale_fill_manual(name = 'Quantiles',
                    values = c('#e5f5e0', '#a1d99b', '#31a354',
                               '#a1d99b', '#e5f5e0')) +
  labs(x = 'Standardised diversity metric ~ fragment size slope estimate',
       y = 'Density') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none',
        legend.justification = c(1, 1),
        text = element_text(size = 7),
        legend.background = element_blank()) 


# this is a version of the plot produced by 7_sim_random_sampling.R
top1 <- cowplot::ggdraw() + 
  cowplot::draw_image(paste0(path2wd, 'intermediate_results/samplingS2000_N40000_mp1.png'),
                                                clip = 'on')

legend_row <- cowplot::plot_grid(legend)


# 3 rows
cowplot::plot_grid(top1, 
                   legend_row,
                   sim_slopes,
                   nrow = 3, 
                   axis = 'tblr',
                   rel_heights = c(1,0.1,2),
                   labels = list('a', '', 'b'),
                   label_size = 8,
                   label_fontface = 'bold')


# set local directory
# figure sized for 2 columns 
# ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision2/figures/Ex_Dat_Fig1.pdf',
#        height = 183, width = 183, units = 'mm')
# 
# ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision3/figures/Ex_Dat_Fig1.png',
#        height = 183, width = 183, units = 'mm')


## calculate the summary stats to report
# slope_coefs %>% 
#   group_by(sigma, metric) %>% 
#   summarise(mean_slope = mean(slope),
#             median_slope = median(slope)) %>% 
#   write.table(file = '~/Dropbox/Frag Database (new)/Manuscript for Nature/revision3/sim_slopes.csv', sep = ',')

# test whether distribution of slopes differ between random and aggregated simulations
dist_test <- slope_coefs %>% 
  select(metric, aggregation, slope) %>% 
  group_by(metric) %>% 
  nest(aggregation, slope) %>% 
  mutate(random_aggr = purrr::map(data, ~t.test(.x %>% 
                                                  filter(aggregation=='Random') %>% 
                                                  select(slope) %>% .$slope,
                                                .x %>% 
                                                  filter(aggregation=='High aggregation') %>% 
                                                  select(slope) %>% .$slope)))
