# get study-level posterior samples from model with only study-level slope variation
# code to plot more detail for each group:
# taxa
library(tidyverse)
library(brms)
library(ggridges)

source('~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/beta_frag_posterior_wrangle.R')


ggplot() +
  facet_grid(continent8 ~ taxa) +
  geom_linerange(data = Jtu_z1i_group_coefs,
                 aes(x = time.since.fragmentation,
                     ymin = Slope_lower, ymax = Slope_upper, 
                     colour = Matrix.category,
                     linetype = time.since.fragmentation,
                     group = interaction(biome, Matrix.category, time.since.fragmentation,dataset_label)),
                 position = position_dodge(width = 1)) +
  geom_point(data = Jtu_z1i_group_coefs,
             aes(x = time.since.fragmentation, y = Slope, colour = Matrix.category,
                 group = interaction(biome, Matrix.category, time.since.fragmentation, dataset_label),
                 shape = biome),
             position = position_dodge(width = 1),
             size = 1.75) +
  geom_hline(yintercept = 0, size = 0.5, lty = 2) +
  geom_hline(data = as.data.frame(Jtu_z1i_fixef),
             aes(yintercept = plogis(Estimate[2]))) +
  geom_rect(data = as.data.frame(Jtu_z1i_fixef),
            aes(xmin = -Inf, xmax = Inf,
                ymin = plogis(Q2.5[2]), ymax = plogis(Q97.5[2])),
            alpha = 0.3) +
  labs(subtitle = expression(paste(beta[Jtu]))) +
  coord_flip() +
  theme_bw()

ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/jtu_study.png',
       width = 290,
       height = 320, units = 'mm')
