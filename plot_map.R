# code to plot map showing locations of studies used in fragmentation analyses
library(tidyverse)

meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

meta <- meta %>% 
  separate(coordinates, into = c('y', 'x'), sep = ', ', remove = F) %>% 
  mutate(x = as.numeric(x),
         y = as.numeric(y))

meta %>% distinct(taxa)

world <- map_data('world') %>% 
  as_tibble()

ggplot() +
  geom_polygon(data=world, 
               aes(long, lat, group = group), colour=NA, fill='#CCCCCC', size=0) +
  geom_point(data = meta,
             aes(x = x, y = y, shape = taxa, colour = taxa),
             position = position_jitter(width = .5)) +
  coord_map('mollweide', ylim = c(-60, 90), xlim = c(-180, 180)) +
  scale_x_continuous(breaks = seq(-180, 180, by = 30)) +
  scale_y_continuous(breaks = c(0, -23.5, 23.5)) +
  scale_shape_manual(name = 'Taxa',
                     values = c('amphibians & reptiles' = 0,
                                'birds' = 1,
                                'plants' = 2,
                                'invertebrates' = 3,
                                'mammals' = 4)) +
  scale_color_manual(name = 'Taxa',
                     values = c('amphibians & reptiles' = '#e41a1c',
                                'birds' = '#984ea3',
                                'plants' = '#4daf4a',
                                'invertebrates' = '#377eb8',
                                'mammals' = '#ff7f00')) +
  labs(x = '', 
       y = '') +
  theme_bw() +
  theme(#panel.grid.major.x = element_blank(), 
        panel.border = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        legend.position = 'top',
        legend.direction = 'horizontal')

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/map.pdf', width = 250, height = 220, units = 'mm')
