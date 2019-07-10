# code to plot results from the models with interactions....
library(tidyverse)
library(brms)
##----load model fits-----
load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_interactions_ref.Rdata')

Sstd2_lognorm_fragSize <- add_criterion(Sstd2_lognorm_fragSize, criterion = 'waic')
Sstd2_ln_fS_tsf <- add_criterion(Sstd2_ln_fS_tsf, criterion = 'waic')
Sstd2_ln_fS_biome <- add_criterion(Sstd2_ln_fS_biome, criterion = 'waic')
Sstd2_ln_fS_matrix <- add_criterion(Sstd2_ln_fS_matrix, criterion = 'waic')
Sstd2_ln_fS_taxa <- add_criterion(Sstd2_ln_fS_taxa, criterion = 'waic')

loo_compare(Sstd2_lognorm_fragSize, Sstd2_ln_fS_tsf, 
            Sstd2_ln_fS_biome, Sstd2_ln_fS_matrix,
            Sstd2_ln_fS_taxa, criterion = 'waic')

# load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_biome_wetland.Rdata')
# for plotting fixed effects
Sstd_fS_matrix_fitted <- cbind(Sstd2_ln_fS_matrix$data,
                             fitted(Sstd2_ln_fS_matrix, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, Matrix.category, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'Matrix.category'))

Sn_fS_matrix_fitted <- cbind(Sn_ln_fS_matrix$data,
                               fitted(Sn_ln_fS_matrix, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, Matrix.category, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'Matrix.category'))

Scov_fS_matrix_fitted <- cbind(Scov_ln_fS_matrix$data,
                             fitted(Scov_ln_fS_matrix, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, Matrix.category, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'Matrix.category'))

Schao_fS_matrix_fitted <- cbind(Schao_ln_fS_matrix$data,
                               fitted(Schao_ln_fS_matrix, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, Matrix.category, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'Matrix.category'))


S_PIE_fS_matrix_fitted <- cbind(S_PIE_ln_fS_matrix$data,
                                fitted(S_PIE_ln_fS_matrix, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, Matrix.category, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'Matrix.category'))

N_std_fS_matrix_fitted <- cbind(N_std_ln_fS_matrix$data,
                                fitted(N_std_ln_fS_matrix, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, Matrix.category, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'Matrix.category'))

# repeat for taxa
Sstd_fS_taxa_fitted <- cbind(Sstd2_ln_fS_taxa$data,
                               fitted(Sstd2_ln_fS_taxa, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, taxa, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'taxa'))

Sn_fS_taxa_fitted <- cbind(Sn_ln_fS_taxa$data,
                             fitted(Sn_ln_fS_taxa, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, taxa, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'taxa'))

Scov_fS_taxa_fitted <- cbind(Scov_ln_fS_taxa$data,
                               fitted(Scov_ln_fS_taxa, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, taxa, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'taxa'))

Schao_fS_taxa_fitted <- cbind(Schao_ln_fS_taxa$data,
                                fitted(Schao_ln_fS_taxa, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, taxa, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'taxa'))


S_PIE_fS_taxa_fitted <- cbind(S_PIE_ln_fS_taxa$data,
                                fitted(S_PIE_ln_fS_taxa, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, taxa, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'taxa'))

N_std_fS_taxa_fitted <- cbind(N_std_ln_fS_taxa$data,
                                fitted(N_std_ln_fS_taxa, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, taxa, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'taxa'))

# repeat for time since fragmentation
Sstd_fS_tsf_fitted <- cbind(Sstd2_ln_fS_tsf$data,
                             fitted(Sstd2_ln_fS_tsf, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, time.since.fragmentation, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'time.since.fragmentation'))

Sn_fS_tsf_fitted <- cbind(Sn_ln_fS_tsf$data,
                           fitted(Sn_ln_fS_tsf, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, time.since.fragmentation, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'time.since.fragmentation'))

Scov_fS_tsf_fitted <- cbind(Scov_ln_fS_tsf$data,
                             fitted(Scov_ln_fS_tsf, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, time.since.fragmentation, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'time.since.fragmentation'))

Schao_fS_tsf_fitted <- cbind(Schao_ln_fS_tsf$data,
                              fitted(Schao_ln_fS_tsf, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, time.since.fragmentation, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'time.since.fragmentation'))


S_PIE_fS_tsf_fitted <- cbind(S_PIE_ln_fS_tsf$data,
                              fitted(S_PIE_ln_fS_tsf, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, time.since.fragmentation, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'time.since.fragmentation'))

N_std_fS_tsf_fitted <- cbind(N_std_ln_fS_tsf$data,
                              fitted(N_std_ln_fS_tsf, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, time.since.fragmentation, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'time.since.fragmentation'))

# repeat for biome
Sstd_fS_biome_fitted <- cbind(Sstd2_ln_fS_biome$data,
                            fitted(Sstd2_ln_fS_biome, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, biome, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'biome'))

Sn_fS_biome_fitted <- cbind(Sn_ln_fS_biome$data,
                          fitted(Sn_ln_fS_biome, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, biome, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'biome'))

Scov_fS_biome_fitted <- cbind(Scov_ln_fS_biome$data,
                            fitted(Scov_ln_fS_biome, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, biome, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'biome'))

Schao_fS_biome_fitted <- cbind(Schao_ln_fS_biome$data,
                             fitted(Schao_ln_fS_biome, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, biome, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'biome'))


S_PIE_fS_biome_fitted <- cbind(S_PIE_ln_fS_biome$data,
                             fitted(S_PIE_ln_fS_biome, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, biome, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'biome'))

N_std_fS_biome_fitted <- cbind(N_std_ln_fS_biome$data,
                             fitted(N_std_ln_fS_biome, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, biome, frag_size_num),
             by = c('dataset_label', 'c.lfs', 'biome'))

## now build the plots: matrix first--------------
Sstd_matrix_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std_2, colour = Matrix.category)) +
  geom_line(data = Sstd_fS_matrix_fitted,
            aes(x = frag_size_num, y = Estimate, colour = Matrix.category),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Sstd_fS_matrix_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = Matrix.category, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[std]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

# get a single legend for all panels
source('~/Dropbox/1current/R_random/functions/gg_legend.R')
matrix_legend <- gg_legend(Sstd_matrix_interaction)

Sn_matrix_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n, colour = Matrix.category)) +
  geom_line(data = Sn_fS_matrix_fitted,
            aes(x = frag_size_num, y = Estimate, colour = Matrix.category),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Sn_fS_matrix_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = Matrix.category, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[n]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

Scov_matrix_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_cov, colour = Matrix.category)) +
  geom_line(data = Scov_fS_matrix_fitted,
            aes(x = frag_size_num, y = Estimate, colour = Matrix.category),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Scov_fS_matrix_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = Matrix.category, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[cov]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

Schao_matrix_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_chao, colour = Matrix.category)) +
  geom_line(data = Schao_fS_matrix_fitted,
            aes(x = frag_size_num, y = Estimate, colour = Matrix.category),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Schao_fS_matrix_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = Matrix.category, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(2, 16, 128, 1024)) +#
  labs(x = '',
       y = expression(paste(S[chao]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

Spie_matrix_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_PIE, colour = Matrix.category)) +
  geom_line(data = S_PIE_fS_matrix_fitted,
            aes(x = frag_size_num, y = Estimate, colour = Matrix.category),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = S_PIE_fS_matrix_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = Matrix.category, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[PIE]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

N_std_matrix_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = N_std, colour = Matrix.category)) +
  geom_line(data = N_std_fS_matrix_fitted,
            aes(x = frag_size_num, y = Estimate, colour = Matrix.category),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = N_std_fS_matrix_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = Matrix.category, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000)) +
  labs(x = '',
       y = expression(paste(N[std]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

matrix_regressions <- cowplot::plot_grid(Sstd_matrix_interaction,
                                         Sn_matrix_interaction,
                                         Scov_matrix_interaction,
                                         Schao_matrix_interaction,
                                         Spie_matrix_interaction,
                                         N_std_matrix_interaction)

cowplot::plot_grid(matrix_legend,
                   matrix_regressions,
                   ncol = 1,
                   rel_heights = c(0.05,1)) +
  cowplot::draw_label("Fragment size (hectares)", x = 0.5, y = 0.01)

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/fS_matrix_interactions.pdf',
#        width = 290, height = 200, units = 'mm')



## now build the plots: taxa--------------
Sstd_taxa_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std_2, colour = taxa)) +
  geom_line(data = Sstd_fS_taxa_fitted,
            aes(x = frag_size_num, y = Estimate, colour = taxa),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Sstd_fS_taxa_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = taxa, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[std]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

# get a single legend for all panels
matrix_legend <- gg_legend(Sstd_taxa_interaction)

Sn_taxa_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n, colour = taxa)) +
  geom_line(data = Sn_fS_taxa_fitted,
            aes(x = frag_size_num, y = Estimate, colour = taxa),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Sn_fS_taxa_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = taxa, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[n]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

Scov_taxa_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_cov, colour = taxa)) +
  geom_line(data = Scov_fS_taxa_fitted,
            aes(x = frag_size_num, y = Estimate, colour = taxa),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Scov_fS_taxa_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = taxa, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[cov]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

Schao_taxa_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_chao, colour = taxa)) +
  geom_line(data = Schao_fS_taxa_fitted,
            aes(x = frag_size_num, y = Estimate, colour = taxa),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Schao_fS_taxa_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = taxa, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(2, 16, 128, 1024)) +#
  labs(x = '',
       y = expression(paste(S[chao]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

Spie_taxa_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_PIE, colour = taxa)) +
  geom_line(data = S_PIE_fS_taxa_fitted,
            aes(x = frag_size_num, y = Estimate, colour = taxa),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = S_PIE_fS_taxa_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = taxa, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[PIE]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

N_std_taxa_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = N_std, colour = taxa)) +
  geom_line(data = N_std_fS_taxa_fitted,
            aes(x = frag_size_num, y = Estimate, colour = taxa),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = N_std_fS_taxa_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = taxa, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000)) +
  labs(x = '',
       y = expression(paste(N[std]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

matrix_regressions <- cowplot::plot_grid(Sstd_taxa_interaction,
                                         Sn_taxa_interaction,
                                         Scov_taxa_interaction,
                                         Schao_taxa_interaction,
                                         Spie_taxa_interaction,
                                         N_std_taxa_interaction)

cowplot::plot_grid(matrix_legend,
                   matrix_regressions,
                   ncol = 1,
                   rel_heights = c(0.05,1)) +
  cowplot::draw_label("Fragment size (hectares)", x = 0.5, y = 0.01)

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/fS_taxa_interactions.pdf',
#        width = 290, height = 200, units = 'mm')



## now build the plots: time--------------
Sstd_tsf_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std_2, colour = time.since.fragmentation)) +
  geom_line(data = Sstd_fS_tsf_fitted,
            aes(x = frag_size_num, y = Estimate, colour = time.since.fragmentation),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Sstd_fS_tsf_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = time.since.fragmentation, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[std]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

# get a single legend for all panels
matrix_legend <- gg_legend(Sstd_tsf_interaction)

Sn_tsf_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n, colour = time.since.fragmentation)) +
  geom_line(data = Sn_fS_tsf_fitted,
            aes(x = frag_size_num, y = Estimate, colour = time.since.fragmentation),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Sn_fS_tsf_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = time.since.fragmentation, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[n]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

Scov_tsf_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_cov, colour = time.since.fragmentation)) +
  geom_line(data = Scov_fS_tsf_fitted,
            aes(x = frag_size_num, y = Estimate, colour = time.since.fragmentation),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Scov_fS_tsf_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = time.since.fragmentation, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[cov]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

Schao_tsf_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_chao, colour = time.since.fragmentation)) +
  geom_line(data = Schao_fS_tsf_fitted,
            aes(x = frag_size_num, y = Estimate, colour = time.since.fragmentation),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Schao_fS_tsf_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = time.since.fragmentation, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(2, 16, 128, 1024)) +#
  labs(x = '',
       y = expression(paste(S[chao]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

Spie_tsf_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_PIE, colour = time.since.fragmentation)) +
  geom_line(data = S_PIE_fS_tsf_fitted,
            aes(x = frag_size_num, y = Estimate, colour = time.since.fragmentation),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = S_PIE_fS_tsf_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = time.since.fragmentation, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[PIE]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

N_std_tsf_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = N_std, colour = time.since.fragmentation)) +
  geom_line(data = N_std_fS_tsf_fitted,
            aes(x = frag_size_num, y = Estimate, colour = time.since.fragmentation),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = N_std_fS_tsf_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = time.since.fragmentation, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000)) +
  labs(x = '',
       y = expression(paste(N[std]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

matrix_regressions <- cowplot::plot_grid(Sstd_tsf_interaction,
                                         Sn_tsf_interaction,
                                         Scov_tsf_interaction,
                                         Schao_tsf_interaction,
                                         Spie_tsf_interaction,
                                         N_std_tsf_interaction)

cowplot::plot_grid(matrix_legend,
                   matrix_regressions,
                   ncol = 1,
                   rel_heights = c(0.05,1)) +
  cowplot::draw_label("Fragment size (hectares)", x = 0.5, y = 0.01)

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/fS_tsf_interactions.pdf',
#        width = 290, height = 200, units = 'mm')



## now build the plots: biome--------------
Sstd_biome_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std_2, colour = biome)) +
  geom_line(data = Sstd_fS_biome_fitted,
            aes(x = frag_size_num, y = Estimate, colour = biome),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Sstd_fS_biome_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = biome, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[std]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

# get a single legend for all panels
matrix_legend <- gg_legend(Sstd_biome_interaction)

Sn_biome_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n, colour = biome)) +
  geom_line(data = Sn_fS_biome_fitted,
            aes(x = frag_size_num, y = Estimate, colour = biome),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Sn_fS_biome_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = biome, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[n]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

Scov_biome_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_cov, colour = biome)) +
  geom_line(data = Scov_fS_biome_fitted,
            aes(x = frag_size_num, y = Estimate, colour = biome),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Scov_fS_biome_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = biome, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[cov]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

Schao_biome_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_chao, colour = biome)) +
  geom_line(data = Schao_fS_biome_fitted,
            aes(x = frag_size_num, y = Estimate, colour = biome),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Schao_fS_biome_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = biome, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(2, 16, 128, 1024)) +#
  labs(x = '',
       y = expression(paste(S[chao]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

Spie_biome_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_PIE, colour = biome)) +
  geom_line(data = S_PIE_fS_biome_fitted,
            aes(x = frag_size_num, y = Estimate, colour = biome),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = S_PIE_fS_biome_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = biome, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  labs(x = '',
       y = expression(paste(S[PIE]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

N_std_biome_interaction <- ggplot() +
  geom_point(data = frag,
             aes(x = frag_size_num, y = N_std, colour = biome)) +
  geom_line(data = N_std_fS_biome_fitted,
            aes(x = frag_size_num, y = Estimate, colour = biome),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = N_std_fS_biome_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5, fill = biome, linetype = NA),
              alpha = 0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000)) +
  labs(x = '',
       y = expression(paste(N[std]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        legend.direction = 'horizontal')

matrix_regressions <- cowplot::plot_grid(Sstd_biome_interaction,
                                         Sn_biome_interaction,
                                         Scov_biome_interaction,
                                         Schao_biome_interaction,
                                         Spie_biome_interaction,
                                         N_std_biome_interaction)

cowplot::plot_grid(matrix_legend,
                   matrix_regressions,
                   ncol = 1,
                   rel_heights = c(0.05,1)) +
  cowplot::draw_label("Fragment size (hectares)", x = 0.5, y = 0.01)

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/fS_biome_interactions.pdf',
#        width = 290, height = 200, units = 'mm')
