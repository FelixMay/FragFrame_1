# code to plot map showing locations of studies used in fragmentation analyses
library(tidyverse)

meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

meta <- meta %>% 
  separate(coordinates, into = c('y', 'x'), sep = ', ', remove = F) %>% 
  mutate(x = as.numeric(x),
         y = as.numeric(y),
         latitude = climate)

meta %>% distinct(taxa)

world <- map_data('world') %>% 
  as_tibble()

map_lat <- ggplot() +
  geom_polygon(data=world, 
               aes(long, lat, group = group), colour=NA, fill='#CCCCCC', size=0) +
  geom_point(data = meta,
             aes(x = x, y = y, shape = taxa, colour = latitude),
             position = position_jitter(width = .5),
             size = 2, stroke = 1.15) +
  coord_map('mollweide', ylim = c(-60, 90), xlim = c(-180, 180)) +
  scale_x_continuous(breaks = seq(-180, 180, by = 30)) +
  scale_y_continuous(breaks = c(0, -23.5, 23.5)) +
  scale_shape_manual(name = '',
                     values = c('amphibians & reptiles' = 0,
                                'birds' = 1,
                                'plants' = 2,
                                'invertebrates' = 3,
                                'mammals' = 4)) +
  scale_color_brewer(name = '', 
                     type = 'qual', palette = 'Set1') +
  labs(x = '', 
       y = '') +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.border = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        legend.position = 'top',
        legend.direction = 'horizontal') +
  guides(shape = guide_legend(nrow = 2))

map_cont7 <- ggplot() +
  geom_polygon(data=world, 
               aes(long, lat, group = group), colour=NA, fill='#CCCCCC', size=0) +
  geom_point(data = meta,
             aes(x = x, y = y, shape = taxa, colour = continent7),
             position = position_jitter(width = .5),
             size = 2, stroke = 1.15) +
  coord_map('mollweide', ylim = c(-60, 90), xlim = c(-180, 180)) +
  scale_x_continuous(breaks = seq(-180, 180, by = 30)) +
  scale_y_continuous(breaks = c(0, -23.5, 23.5)) +
  scale_shape_manual(guide = F, 
                     name = 'Taxa',
                     values = c('amphibians & reptiles' = 0,
                                'birds' = 1,
                                'plants' = 2,
                                'invertebrates' = 3,
                                'mammals' = 4)) +
  scale_color_brewer(name = '', 
                     type = 'qual', palette = 'Set2') +
  labs(x = '', 
       y = '') +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.border = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        legend.position = 'top',
        legend.direction = 'horizontal')

map_cont8 <- ggplot() +
  geom_polygon(data=world, 
               aes(long, lat, group = group), colour=NA, fill='#CCCCCC', size=0) +
  geom_point(data = meta,
             aes(x = x, y = y, shape = taxa, colour = continent8),
             position = position_jitter(width = .5),
             size = 2, stroke = 1.15) +
  coord_map('mollweide', ylim = c(-60, 90), xlim = c(-180, 180)) +
  scale_x_continuous(breaks = seq(-180, 180, by = 30)) +
  scale_y_continuous(breaks = c(0, -23.5, 23.5)) +
  scale_shape_manual(guide = F, 
                     name = 'Taxa',
                     values = c('amphibians & reptiles' = 0,
                                'birds' = 1,
                                'plants' = 2,
                                'invertebrates' = 3,
                                'mammals' = 4)) +
  scale_color_brewer(name = '', 
                     type = 'qual', palette = 'Set3') +
  labs(x = '', 
       y = '') +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.border = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        legend.position = 'top',
        legend.direction = 'horizontal')

cowplot::plot_grid(map_lat, map_cont7, map_cont8, NULL, nrow = 2)
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/alt_maps.pdf', width = 290, height = 220, units = 'mm')
