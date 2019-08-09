# code to plot map showing locations of studies used in fragmentation analyses

# need to run 0_init_dirs_load_packages.R to load packages and setwd etc

# load the meta data
meta <- read.csv(paste0(path2meta, 'new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

meta <- meta %>% 
  separate(coordinates, into = c('y', 'x'), sep = ', ', remove = F) %>% 
  mutate(x = as.numeric(x),
         y = as.numeric(y),
         latitude = climate)

meta %>% distinct(taxa)
# fix some labels for plotting 
meta$taxa <- factor(meta$taxa,
                         levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                         labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

world <- map_data('world') %>% 
  as_tibble()

map_taxa <- ggplot() +
  geom_polygon(data=world, 
               aes(long, lat, group = group), colour=NA, fill='#CCCCCC', size=0) +
  geom_point(data = meta,
             aes(x = x, y = y, colour = taxa),#,shape = taxa),
             # position = position_jitter(width = .5),
             alpha = 0.6,
             size = 1.5) +
  coord_map('mollweide', ylim = c(-60, 90), xlim = c(-180, 180)) +
  scale_x_continuous(breaks = seq(-180, 180, by = 30)) +
  scale_y_continuous(breaks = c(0, -23.5, 23.5, -60, 60)) +
  # scale_shape_manual(name = 'Taxon group',
  #                    values = c('Amphibians & reptiles' = 0,
  #                               'Birds' = 1,
  #                               'Plants' = 2,
  #                               'Invertebrates' = 3,
  #                               'Mammals' = 4),
  #                    labels = c('Amphibians & reptiles', 'Birds',
  #                               'Plants', 'Invertebrates', 'Mammals')) +
  scale_color_brewer(name = 'Taxon group', 
                     type = 'qual', palette = 'Dark2',
                     labels = c('Amphibians & reptiles', 'Birds',
                                'Plants', 'Invertebrates', 'Mammals')) +
  labs(x = '', 
       y = '') +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = 'black', size = 0.1), 
        panel.border = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        legend.position = 'top',
        legend.direction = 'horizontal',
        plot.margin = unit(c(0,0,0,0), units = 'mm')) +
  guides(shape = guide_legend(nrow = 2),
         colour = guide_legend(nrow = 2)
         )


# figure 1
# we want the conceptual figure on the top panel
top <- cowplot::ggdraw() + cowplot::draw_image('~/Dropbox/Frag Database (new)/analysis_apr19/figures/fig1_top.png',
                                               clip = 'on')
cowplot::plot_grid(top, map_taxa,
                   nrow = 2, align = 'hv',
                   labels = 'auto',
                   rel_heights = c(0.5,1),
                   rel_widths = c(0.5,1))

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/fig1_solidPoints.pdf',
#        width = 150, height = 150, units = 'mm')

