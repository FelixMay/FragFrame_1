# visual inspection of data for fragmentation synthesis
library(tidyverse)

frag <- read_csv('~/Dropbox/Habitat loss meta-analysis/analysis/diversity_metadata.csv')

# we are mostly interested in diversity as a function of entity size
ggplot() +
  geom_point(data = frag,
             aes(x = entity.size, y = S_n)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(2,32,64,128, 256)) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = entity.size, y = S_n),
              colour = 'black') +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[n]))) +
  theme_bw()

# what does the study-level variation look like?
ggplot() +
  geom_point(data = frag,
             aes(x = entity.size, y = S_n)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(2,32,64,128, 256)) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = entity.size, y = S_n),
              colour = 'black') +
  stat_smooth(data = frag,
              method = 'lm', se = F,
              aes(x = entity.size, y = S_n,
                  group = filename, colour = taxa),
              lwd = 0.3) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[n]))) +
  theme_bw()

# How do we want to examine taxa? Within studies, or across all studies.
# are there really studies where we don't know the taxa? check NAs here...
# this plot is taxa across all studies
ggplot() +
  geom_point(data = frag,
             aes(x = entity.size, y = S_n, 
                 colour = filename)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(2,32,64,128, 256)) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = entity.size, y = S_n)) +
  stat_smooth(data = frag,
              method = 'lm', se = F,
              aes(x = entity.size, y = S_n, 
                  colour = taxa)) +
  theme(legend.position = 'none')

# taxa may respond differently depending on the vegetation within the fragments
ggplot() +
  facet_wrap(~veg.fragment) +
  geom_point(data = frag,
             aes(x = entity.size, y = S_n, 
                 colour = taxa)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(2,32,64,128, 256)) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = entity.size, y = S_n)) +
  stat_smooth(data = frag,
              method = 'lm', se = F,
              aes(x = entity.size, y = S_n, 
                  colour = taxa)) +
  theme_bw()

# what about the character of the matrix?
# maybe not enough data for this: smaller range of entity sizes in 'light' matrix
ggplot() +
  facet_wrap(~matrix.category) +
  geom_point(data = frag,
             aes(x = entity.size, y = S_n, 
                 colour = taxa)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(2,32,64,128, 256)) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = entity.size, y = S_n)) +
  stat_smooth(data = frag,
              method = 'lm', se = F,
              aes(x = entity.size, y = S_n, 
                  colour = taxa)) +
  theme_bw()

# time since fragmentation
# long time since fragmentation has an overall negative slope...
ggplot() +
  facet_wrap(~time.since.fragmentation) +
  geom_point(data = frag,
             aes(x = entity.size, y = S_n,
                 colour = taxa)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log2', breaks = c(2,32,64,128, 256)) +
  stat_smooth(data = frag,
              method = 'lm',
              aes(x = entity.size, y = S_n)) +
  stat_smooth(data = frag,
              method = 'lm', se = F,
              aes(x = entity.size, y = S_n,
                  colour = taxa)) +
  theme_bw()
