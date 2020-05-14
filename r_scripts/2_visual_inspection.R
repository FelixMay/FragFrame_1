# visual inspection of data for fragmentation synthesis

# need to run 0_init_dirs_load_packages.R for packages and directory paths

frag <- read_csv(paste0(path2wd, 'intermediate_results/2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# load the meta data
meta <- read_delim(paste0(path2wd, 'data/new_meta_2_merge.csv'), delim = ';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

# and join 
frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')


# we are mostly interested in diversity as a function of fragment size
ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_PIE_mean, colour = 'S_PIE'),
             alpha = 0.5) +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n_mean, colour = 'S_n'),
             alpha = 0.5) +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_cov_mean, colour = 'S_cov'),
             alpha = 0.5) +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std_mean, colour = 'S_std_2'),
             alpha = 0.5) +
  stat_smooth(data = frag,
              method = 'lm', 
              aes(x = frag_size_num, y = S_PIE_mean,
                  colour = 'S_PIE')) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = frag_size_num, y = S_n_mean, colour = 'S_n')) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = frag_size_num, y = S_cov_mean, colour = 'S_cov')) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = frag_size_num, y = S_std_mean, colour = 'S_std_2')) +
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


# what does the study-level variation look like?
ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n_mean, colour = sample_design)) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = frag_size_num, y = S_n_mean),
              colour = 'black') +
  stat_smooth(data = frag,
              method = 'lm', se = F,
              aes(x = frag_size_num, y = S_n_mean,
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

# one example of another covariate of interest: taxa
ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n_mean, 
                 colour = taxa),
             alpha = 0.3) +
  stat_smooth(data = frag,
              method = 'lm', se = F,
              aes(x = frag_size_num, y = S_n_mean, 
                  colour = taxa)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(2,32,64,128,256)) +
  scale_colour_viridis_d(name = 'Taxa') +
  theme_bw() +
  theme(legend.position = c(0.1, 0.9))

