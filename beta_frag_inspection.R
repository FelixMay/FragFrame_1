# code to inspect beta-diversity for fragmentation synthesis (two-scales in data: fragment and study)
# code to run on rstudio server
rm(list=ls())
library(tidyverse)

# load the data
frag_beta <- read_csv('~/Dropbox/Frag Database (new)/files_datapaper/Analysis/2_betapart_frag_fcont_10_mabund_as_is.csv')
study_beta <- read_csv('~/Dropbox/Frag Database (new)/files_datapaper/Analysis/2_betapart_study_fcont_10_mabund_as_is.csv')
meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

str(frag_beta)

frag_beta %>% distinct(method)
frag_beta <- frag_beta %>% 
  separate(method, c('family', 'metric'), sep = ', ', remove= F) %>% 
  # add a dirty fix for the zeroes and ones (to facilitate logit-transform)
  # how to model these data? 
  mutate(rch_01 = ifelse(rich==0, rich+1e-3,
                          ifelse(rich==1, rich-1e-6, rich)),
         rpl_01 = ifelse(repl==0, repl+1e-3,
                          ifelse(repl==1, repl-1e-6, repl)))


# compare baselga and podani: what a mess! 
ggplot() +
  facet_grid(metric~family) +
  # plot the data: richness (nestedness) component
  geom_point(data = frag_beta %>% filter(metric!='Sorensen' & metric!='percentage difference'),
                 aes(frag_size_num.x/frag_size_num.y, y = rich,
                     colour = 'rich'),
                 size = 0.1, alpha = 0.1) +
  # stat_density2d(data = frag_beta %>% dplyr::filter(metric!='Sorensen' & metric!='percentage difference'),
  #                aes(frag_size_num.x/frag_size_num.y, y = rich,
  #                    linetype = metric, colour = 'rich'),
  #                alpha = 0.7) +
  # plot the dataL: turnover (replacement) component
  geom_point(data = frag_beta %>% dplyr::filter(metric!='Sorensen' & metric!='percentage difference'),
             aes(frag_size_num.x/frag_size_num.y, y = repl,
                 colour = 'repl'),
             size = 0.1, alpha = 0.1) +
  # stat_density2d(data = frag_beta %>% dplyr::filter(metric!='Sorensen' & metric!='percentage difference'),
  #                aes(frag_size_num.x/frag_size_num.y, y = repl,
  #                    # linetype = metric, 
  #                    colour = 'repl'),
  #                alpha = 0.7) +
  # linear model of Jaccard nestedness
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Jaccard'),
              aes(x = frag_size_num.x/frag_size_num.y, y = rich,
                  # linetype = metric, 
                  colour = 'rich'),
              method = 'lm', se = F) +
  # linear model of Jaccard nestedness: study-level 
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Jaccard'),
              aes(x = frag_size_num.x/frag_size_num.y, y = rich, group = dataset_label,
                  # linetype = metric, 
                  colour = 'rich'),
              size = 0.1, alpha = 0.1,
              method = 'lm', se = F) +
  # linear model of Baselga's Ruzicka metric: nestedness
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Ruzicka'),
              aes(x = frag_size_num.x/frag_size_num.y, y = rich,
                  # linetype = metric, 
                  colour = 'rich'),
              method = 'lm', se = F) +
  # linear model of Baselga's Ruzicka metric: nestedness study-level
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Ruzicka'),
              aes(x = frag_size_num.x/frag_size_num.y, y = rich, group = dataset_label,
                  # linetype = metric, 
                  colour = 'rich'),
              size = 0.1, alpha = 0.1,
              method = 'lm', se = F) +
  # linear model of Baselga's Jaccard turnover metric
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Jaccard'),
              aes(x = frag_size_num.x/frag_size_num.y, y = repl,
                  # linetype = metric, 
                  colour = 'repl'),
              method = 'lm', se = F) +
  # now study-level
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Jaccard'),
              aes(x = frag_size_num.x/frag_size_num.y, y = repl, group = dataset_label,
                  # linetype = metric, 
                  colour = 'repl'),
              size = 0.1, alpha = 0.1,
              method = 'lm', se = F) +
  # linear model of Baselga's Ruzicka turnover metric
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Ruzicka'),
              aes(x = frag_size_num.x/frag_size_num.y, y = repl,
                  # linetype = metric, 
                  colour = 'repl'),
              method = 'lm', se = F) +
  # linear model of Baselga's Ruzicka turnover metric, study-level
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Ruzicka'),
              aes(x = frag_size_num.x/frag_size_num.y, y = repl, group = dataset_label,
                  # linetype = metric, 
                  colour = 'repl'),
              size = 0.1, alpha = 0.1,
              method = 'lm', se = F) +
  # stat_smooth(data = frag_beta %>% filter(method=='Baselga family, percentage difference'),
  #             aes(x = frag_size_num.x/frag_size_num.y, y = rich, 
  #                 linetype = metric, colour = 'rich'),
  #             method = 'lm', se = F) +
  # stat_smooth(data = frag_beta %>% filter(method=='Baselga family, percentage difference'),
  #             aes(x = frag_size_num.x/frag_size_num.y, y = repl, 
  #                 linetype = metric, colour = 'repl'),
  #             method = 'lm', se = F) +
  # # add the podani 'equivalents'
  stat_smooth(data = frag_beta %>% filter(method=='Podani family, Jaccard'),
              aes(x = frag_size_num.x/frag_size_num.y, y = rich,
                  # linetype = metric, 
                  colour = 'rich'),
              method = 'lm', se = F) +
  # # add the podani 'equivalents'
  stat_smooth(data = frag_beta %>% filter(method=='Podani family, Jaccard'),
              aes(x = frag_size_num.x/frag_size_num.y, y = rich, group = dataset_label,
                  # linetype = metric, 
                  colour = 'rich'),
              size = 0.1, alpha = 0.1,
              method = 'lm', se = F) +
  #
  stat_smooth(data = frag_beta %>% filter(method=='Podani family, Ruzicka'),
              aes(x = frag_size_num.x/frag_size_num.y, y = rich,
                  # linetype = metric, 
                  colour = 'rich'),
              method = 'lm', se = F) +
  #
  stat_smooth(data = frag_beta %>% filter(method=='Podani family, Ruzicka'),
              aes(x = frag_size_num.x/frag_size_num.y, y = rich, group = dataset_label,
                  # linetype = metric, 
                  colour = 'rich'),
              size = 0.1, alpha = 0.1,
              method = 'lm', se = F) +
  #
  stat_smooth(data = frag_beta %>% filter(method=='Podani family, Jaccard'),
              aes(x = frag_size_num.x/frag_size_num.y, y = repl,
                  # linetype = metric, 
                  colour = 'repl'),
              method = 'lm', se = F) +
  #
  stat_smooth(data = frag_beta %>% filter(method=='Podani family, Jaccard'),
              aes(x = frag_size_num.x/frag_size_num.y, y = repl, group = dataset_label,
                  # linetype = metric, 
                  colour = 'repl'),
              size = 0.1, alpha = 0.1,
              method = 'lm', se = F) +
  #
  stat_smooth(data = frag_beta %>% filter(method=='Podani family, Ruzicka'),
              aes(x = frag_size_num.x/frag_size_num.y, y = repl,
                  # linetype = metric, 
                  colour = 'repl'),
              method = 'lm', se = F) +
  #
  stat_smooth(data = frag_beta %>% filter(method=='Podani family, Ruzicka'),
              aes(x = frag_size_num.x/frag_size_num.y, y = repl, group = dataset_label,
                  # linetype = metric, 
                  colour = 'repl'),
              size = 0.1, alpha = 0.1,
              method = 'lm', se = F) +
  # stat_smooth(data = frag_beta %>% filter(method=='Podani family, percentage difference'),
  #             aes(x = frag_size_num.x/frag_size_num.y, y = rich, 
  #                 linetype = metric, colour = 'rich'),
  #             method = 'lm', se = F) +
  # stat_smooth(data = frag_beta %>% filter(method=='Podani family, percentage difference'),
  #             aes(x = frag_size_num.x/frag_size_num.y, y = repl,
  #                 linetype = metric, colour = 'repl'),
  #             method = 'lm', se = F) +
  scale_x_continuous(trans =  'log10') +
  # scale_y_continuous(trans = 'logit') +
  scale_color_manual(name = 'dissimilarity component',
                     values = c('repl' = 'blue', 'rich' = 'black')) +
  labs(y = 'Pairwise dissimilarity',
       x = 'Ratio of fragment sizes (log-scale)') +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        panel.grid = element_blank())

ggsave('~/Dropbox/1current/fragmentation_synthesis/results/figures/inspection/beta_frag_family_inspection.png', 
       width = 200, height = 200, units = 'mm')
# baselga and podani are different: following Baselga & Leprieur MEE, I chose to present Baselga
ggplot() +
  # facet_grid(metric~family) +
  # plot the data: richness (nestedness) component
  # geom_point(data = frag_beta %>% filter(family=='Baselga family'),
  #            aes(frag_size_num.x/frag_size_num.y, y = rich,
  #                colour = 'rich'),
  #            size = 0.1, alpha = 0.1) +
  # # stat_density2d(data = frag_beta %>% dplyr::filter(metric!='Sorensen' & metric!='percentage difference'),
  # #                aes(frag_size_num.x/frag_size_num.y, y = rich,
  # #                    linetype = metric, colour = 'rich'),
  # #                alpha = 0.7) +
  # # plot the dataL: turnover (replacement) component
  # geom_point(data = frag_beta %>% filter(family=='Baselga family'),
  #            aes(frag_size_num.x/frag_size_num.y, y = repl,
  #                colour = 'repl'),
  #            size = 0.1, alpha = 0.1) +
  # stat_density2d(data = frag_beta %>% dplyr::filter(metric!='Sorensen' & metric!='percentage difference'),
  #                aes(frag_size_num.x/frag_size_num.y, y = repl,
  #                    # linetype = metric, 
  #                    colour = 'repl'),
  #                alpha = 0.7) +
  # linear model of Jaccard nestedness
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Jaccard'),
              aes(x = frag_size_num.x/frag_size_num.y, y = rich,
                  linetype = metric,
                  colour = 'rich'),
              method = 'lm', se = F) +
  # linear model of Jaccard nestedness: study-level 
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Jaccard'),
              aes(x = frag_size_num.x/frag_size_num.y, y = rich, group = dataset_label,
                  linetype = metric,
                  colour = 'rich'),
              size = 0.1, alpha = 0.1,
              method = 'lm', se = F) +
  # linear model of Baselga's Ruzicka metric: nestedness
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Ruzicka'),
              aes(x = frag_size_num.x/frag_size_num.y, y = rich,
                  linetype = metric,
                  colour = 'rich'),
              method = 'lm', se = F) +
  # linear model of Baselga's Ruzicka metric: nestedness study-level
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Ruzicka'),
              aes(x = frag_size_num.x/frag_size_num.y, y = rich, group = dataset_label,
                  linetype = metric,
                  colour = 'rich'),
              size = 0.1, alpha = 0.1,
              method = 'lm', se = F) +
  # linear model of Baselga's Jaccard turnover metric
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Jaccard'),
              aes(x = frag_size_num.x/frag_size_num.y, y = repl,
                  linetype = metric,
                  colour = 'repl'),
              method = 'lm', se = F) +
  # now study-level
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Jaccard'),
              aes(x = frag_size_num.x/frag_size_num.y, y = repl, group = dataset_label,
                  linetype = metric,
                  colour = 'repl'),
              size = 0.1, alpha = 0.1,
              method = 'lm', se = F) +
  # linear model of Baselga's Ruzicka turnover metric
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Ruzicka'),
              aes(x = frag_size_num.x/frag_size_num.y, y = repl,
                  linetype = metric,
                  colour = 'repl'),
              method = 'lm', se = F) +
  # linear model of Baselga's Ruzicka turnover metric, study-level
  stat_smooth(data = frag_beta %>% filter(method=='Baselga family, Ruzicka'),
              aes(x = frag_size_num.x/frag_size_num.y, y = repl, group = dataset_label,
                  linetype = metric,
                  colour = 'repl'),
              size = 0.1, alpha = 0.1,
              method = 'lm', se = F) +
    scale_x_continuous(trans =  'log10') +
  # scale_y_continuous(trans = 'logit') +
  scale_color_manual(name = 'dissimilarity component',
                     values = c('repl' = 'blue', 'rich' = 'black')) +
  scale_linetype_manual(name = 'metric',
                     values = c('Jaccard' = 1, 'Ruzicka' = 2)) +
  labs(y = 'Pairwise dissimilarity') +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        panel.grid = element_blank())

hist(frag_beta$log10_ratio_area - median(frag_beta$log10_ratio_area))
