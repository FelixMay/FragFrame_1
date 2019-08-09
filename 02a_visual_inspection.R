# visual inspection of data for fragmentation synthesis

# need to run 0_init_dirs_load_packages.R for packages and directory paths

frag <- read_csv(paste0(path2data, '2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# load the meta data
meta <- read.csv(paste0(path2meta, '/new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

# and join 
frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')

# if you want to save these figures...
setwd(paste(path2temp,"figs/visual_inspection/", sep = ""))

# we are mostly interested in diversity as a function of fragment size
ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_PIE, colour = 'S_PIE'),
             alpha = 0.5) +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n, colour = 'S_n'),
             alpha = 0.5) +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_cov, colour = 'S_cov'),
             alpha = 0.5) +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std_2, colour = 'S_std_2'),
             alpha = 0.5) +
  stat_smooth(data = frag,
              method = 'lm', 
              aes(x = frag_size_num, y = S_PIE,
                  colour = 'S_PIE')) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = frag_size_num, y = S_n, colour = 'S_n')) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = frag_size_num, y = S_cov, colour = 'S_cov')) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = frag_size_num, y = S_std_2, colour = 'S_std_2')) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(2,32,64,128, 256)) +
  scale_colour_manual(name = 'metric',
                      values = c('N_std' = '#00455c', 'S_PIE' = '#006a6f', 
                                 'S_n' = '#008c57', 'S_cov' = '#83a51b', 
                                 'S_std_2' = '#ffa600')) +
  labs(x = 'Fragment size (hectares)',
       y = 'Richness (standardised)') +
  theme_bw() +
  theme(legend.key = element_blank(),
        legend.position = c(0.1,0.9)) +
  guides(colour = guide_legend(override.aes = list(fill=NA)))

# ggsave('standardised_S_fragmentSize.pdf', width = 290, height = 200, units = 'mm')

# Total and standardized N
ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = N_std, colour = 'N_std'),
             alpha = 0.5) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = frag_size_num, y = N_std,
                  colour = 'N_std')) +
   scale_x_continuous(trans = 'log10') +
   scale_y_continuous(trans = 'log10') +
  scale_colour_manual(name = 'metric', guide = F,
                      values = c('N_std' = '#00455c')) +
   labs(x = 'Fragment size (hectares)',
        y = 'No. of individuals') +
   theme_bw() 
   
# ggsave('N_fragmentSize.pdf', width = 290, height = 200, units = 'mm')


# what does the study-level variation look like?
ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n, colour = sample_design)) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = frag_size_num, y = S_n),
              colour = 'black') +
  stat_smooth(data = frag,
              method = 'lm', se = F,
              aes(x = frag_size_num, y = S_n,
                  group = dataset_label,
                  colour = sample_design),
              lwd = 0.3) +
  scale_colour_viridis_d(name = 'Sample design') +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(4,32,64,128,256)) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[n]))) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.9))

# ggsave('S_n_fragmentSize_studyLevel.pdf', width = 290, height = 200, units = 'mm')

# How do we want to examine taxa? Within studies, or across all studies?
# this plot is taxa across all studies...
ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n, 
                 colour = taxa),
             alpha = 0.3) +
  stat_smooth(data = frag,
              method = 'lm', se = F,
              aes(x = frag_size_num, y = S_n, 
                  colour = taxa)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(2,32,64,128,256)) +
  scale_colour_viridis_d(name = 'Taxa') +
  theme_bw() +
  theme(legend.position = c(0.1, 0.9))

#....but I think we should allow the taxa to vary within studies too:
ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n, 
                 colour = taxa), 
             alpha = 0.3, size = 1.5) +
  stat_smooth(data = frag,
              method = 'lm', se = F,
              aes(x = frag_size_num, y = S_n,
                  colour = taxa),
              size = 1.5) +
  stat_smooth(data = frag,
              method = 'lm', se = F,
              aes(x = frag_size_num, y = S_n, group = dataset_label,
                  colour = taxa),
              size = 0.3) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(2,32,64,128, 256)) +
  scale_colour_viridis_d(name = 'Taxa') +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[n]))) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.8),
        legend.background = element_blank())

# ggsave('standardised_S_fragmentSize_x_taxa_studyLevel.pdf', width = 290, height = 200, units = 'mm')

# what about the character of the matrix?
# add factor for ordering facets
frag$f.Matrix.category <- factor(frag$Matrix.category, levels = c('light filter', 'intermediate', 'harsh filter')) 

ggplot() +
  facet_wrap(~f.Matrix.category) +
  # throw out the fragments for which we don't have Matrix.category data
  geom_point(data = frag %>% filter(!is.na(Matrix.category)),
             aes(x = frag_size_num, y = S_n, 
                 colour = taxa)) +
  stat_smooth(data = frag %>% filter(!is.na(Matrix.category)),
              method = 'lm',
              aes(x = frag_size_num, y = S_n),
              colour = 'black') +
  stat_smooth(data = frag %>% filter(!is.na(Matrix.category)),
              method = 'lm', se = F,
              aes(x = frag_size_num, y = S_n, group = dataset_label,
                  colour = taxa),
              lwd = 0.3) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(2,32,64,128, 256)) +
  scale_colour_viridis_d(name = 'Taxa') +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[n]))) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.8),
        legend.background = element_blank())

# ggsave('Sn_fragmentSize_x_taxa_x_matrixCategory_studyLevel.pdf', width = 290, height = 140, units = 'mm')

# time since fragmentation
# order for facets...
frag$f.tsf <- factor(frag$time.since.fragmentation, 
                     levels = c('Recent (less than 20 years)', 'Intermediate (20-100 years)', 'long (100+ years)'),
                     labels = c('< 20 years', '20 - 100 years', '> 100 years'))

# long time since fragmentation looks to have a flatter (negative?) slope
ggplot() +
  facet_wrap(~f.tsf) +
  # throw out the fragments for which we don't have Matrix.category data
  geom_point(data = frag %>% filter(!is.na(time.since.fragmentation)),
             aes(x = frag_size_num, y = S_n, 
                 colour = taxa)) +
  stat_smooth(data = frag %>% filter(!is.na(time.since.fragmentation)),
              method = 'lm',
              aes(x = frag_size_num, y = S_n),
              colour = 'black') +
  stat_smooth(data = frag %>% filter(!is.na(time.since.fragmentation)),
              method = 'lm', se = F,
              aes(x = frag_size_num, y = S_n, group = dataset_label,
                  colour = taxa),
              lwd = 0.3) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(2,32,64,128, 256)) +
  scale_colour_viridis_d(name = 'Taxa') +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[n]))) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.9),
        legend.background = element_blank()) +
  guides(colour = guide_legend(ncol = 2))

# ggsave('Sn_fragmentSize_x_taxa_x_timeSinceFrag_studyLevel.pdf', width = 290, height = 140, units = 'mm')