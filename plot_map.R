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
             aes(x = x, y = y)) +
  coord_map('mollweide', ylim = c(-60, 90), xlim = c(-180, 180)) +
  labs(x = '', 
       y = '') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank())
