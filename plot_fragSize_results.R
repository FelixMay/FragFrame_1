# code to plot results of fragment size (only) models
library(tidyverse)
load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms.Rdata')

#------wrangle for plotting
# for plotting fixed effects----------------
Sstd1_fS_fitted <- cbind(Sstd1_lognorm_fragSize$data,
                               fitted(Sstd1_lognorm_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, frag_size_num),
             by = c('dataset_label', 'c.lfs'))

Sstd1_lognorm_fragSize_fixef <- fixef(Sstd1_lognorm_fragSize)

Sstd2_fS_fitted <- cbind(Sstd2_lognorm_fragSize$data,
                         fitted(Sstd2_lognorm_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, frag_size_num),
             by = c('dataset_label', 'c.lfs'))

Sstd2_lognorm_fragSize_fixef <- fixef(Sstd2_lognorm_fragSize)

Sn_fS_fitted <- cbind(Sn_lognorm_fragSize$data,
                             fitted(Sn_lognorm_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, frag_size_num),
             by = c('dataset_label', 'c.lfs'))

Sn_lognorm_fragSize_fixef <- fixef(Sn_lognorm_fragSize)

Scov_fS_fitted <- cbind(Scov_lognorm_fragSize$data,
                      fitted(Scov_lognorm_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, frag_size_num),
             by = c('dataset_label', 'c.lfs'))

Scov_lognorm_fragSize_fixef <- fixef(Scov_lognorm_fragSize)

Schao_fS_fitted <- cbind(S_chao_lognorm_fragSize$data,
                        fitted(S_chao_lognorm_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, frag_size_num),
             by = c('dataset_label', 'c.lfs'))

Schao_lognorm_fragSize_fixef <- fixef(S_chao_lognorm_fragSize)

S_PIE_fS_fitted <- cbind(S_PIE_lognorm_fragSize$data,
                                fitted(S_PIE_lognorm_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, frag_size_num),
             by = c('dataset_label', 'c.lfs'))

S_PIE_lognorm_fragSize_fixef <- fixef(S_PIE_lognorm_fragSize)

Nstd_fS_fitted <- cbind(Nstd_lognorm_fragSize$data,
                               fitted(Nstd_lognorm_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, frag_size_num),
             by = c('dataset_label', 'c.lfs'))

Nstd_lognorm_fragSize_fixef <- fixef(Nstd_lognorm_fragSize)


# for plotting the random-effects----------------
Sstd1_lognorm_fragSize_coef <- coef(Sstd1_lognorm_fragSize)
Sstd2_lognorm_fragSize_coef <- coef(Sstd2_lognorm_fragSize)
Sn_lognorm_fragSize_coef <- coef(Sn_lognorm_fragSize)
Scov_lognorm_fragSize_coef <- coef(Scov_lognorm_fragSize)
Schao_lognorm_fragSize_coef <- coef(S_chao_lognorm_fragSize)
S_PIE_fS_coef <- coef(S_PIE_lognorm_fragSize)
Nstd_fS_coef <- coef(Nstd_lognorm_fragSize)

Sstd1_lognorm_fragSize_group_coefs <- bind_cols(Sstd1_lognorm_fragSize_coef[[1]][,,'Intercept'] %>% 
                                          as_tibble() %>% 
                                          mutate(Intercept = Estimate,
                                                 Intercept_lower = Q2.5,
                                                 Intercept_upper = Q97.5,
                                                 dataset_label = rownames(Sstd1_lognorm_fragSize_coef[[1]][,,'Intercept'])) %>% 
                                          dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                        Sstd1_lognorm_fragSize_coef[[1]][,,'c.lfs'] %>% 
                                          as_tibble() %>% 
                                          mutate(Slope = Estimate,
                                                 Slope_lower = Q2.5,
                                                 Slope_upper = Q97.5) %>% 
                                          dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(dataset_label) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'dataset_label')

Sstd2_lognorm_fragSize_group_coefs <- bind_cols(Sstd2_lognorm_fragSize_coef[[1]][,,'Intercept'] %>% 
                                                  as_tibble() %>% 
                                                  mutate(Intercept = Estimate,
                                                         Intercept_lower = Q2.5,
                                                         Intercept_upper = Q97.5,
                                                         dataset_label = rownames(Sstd2_lognorm_fragSize_coef[[1]][,,'Intercept'])) %>% 
                                                  dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                                Sstd2_lognorm_fragSize_coef[[1]][,,'c.lfs'] %>% 
                                                  as_tibble() %>% 
                                                  mutate(Slope = Estimate,
                                                         Slope_lower = Q2.5,
                                                         Slope_upper = Q97.5) %>% 
                                                  dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(dataset_label) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'dataset_label')

Sn_lognorm_fragSize_group_coefs <- bind_cols(Sn_lognorm_fragSize_coef[[1]][,,'Intercept'] %>% 
                                        as_tibble() %>% 
                                        mutate(Intercept = Estimate,
                                               Intercept_lower = Q2.5,
                                               Intercept_upper = Q97.5,
                                               dataset_label = rownames(Sn_lognorm_fragSize_coef[[1]][,,'Intercept'])) %>% 
                                        dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                      Sn_lognorm_fragSize_coef[[1]][,,'c.lfs'] %>% 
                                        as_tibble() %>% 
                                        mutate(Slope = Estimate,
                                               Slope_lower = Q2.5,
                                               Slope_upper = Q97.5) %>% 
                                        dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(dataset_label) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'dataset_label')

Scov_lognorm_fragSize_group_coefs <- bind_cols(Scov_lognorm_fragSize_coef[[1]][,,'Intercept'] %>% 
                                               as_tibble() %>% 
                                               mutate(Intercept = Estimate,
                                                      Intercept_lower = Q2.5,
                                                      Intercept_upper = Q97.5,
                                                      dataset_label = rownames(Scov_lognorm_fragSize_coef[[1]][,,'Intercept'])) %>% 
                                               dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                               Scov_lognorm_fragSize_coef[[1]][,,'c.lfs'] %>% 
                                               as_tibble() %>% 
                                               mutate(Slope = Estimate,
                                                      Slope_lower = Q2.5,
                                                      Slope_upper = Q97.5) %>% 
                                               dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(dataset_label) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'dataset_label')

Schao_lognorm_fragSize_group_coefs <- bind_cols(Schao_lognorm_fragSize_coef[[1]][,,'Intercept'] %>% 
                                                 as_tibble() %>% 
                                                 mutate(Intercept = Estimate,
                                                        Intercept_lower = Q2.5,
                                                        Intercept_upper = Q97.5,
                                                        dataset_label = rownames(Schao_lognorm_fragSize_coef[[1]][,,'Intercept'])) %>% 
                                                 dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                                Schao_lognorm_fragSize_coef[[1]][,,'c.lfs'] %>% 
                                                 as_tibble() %>% 
                                                 mutate(Slope = Estimate,
                                                        Slope_lower = Q2.5,
                                                        Slope_upper = Q97.5) %>% 
                                                 dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(dataset_label) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'dataset_label')

S_PIE_fragSize_group_coefs <- bind_cols(S_PIE_fS_coef[[1]][,,'Intercept'] %>% 
                                           as_tibble() %>% 
                                           mutate(Intercept = Estimate,
                                                  Intercept_lower = Q2.5,
                                                  Intercept_upper = Q97.5,
                                                  dataset_label = rownames(S_PIE_fS_coef[[1]][,,'Intercept'])) %>% 
                                           dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                        S_PIE_fS_coef[[1]][,,'c.lfs'] %>% 
                                           as_tibble() %>% 
                                           mutate(Slope = Estimate,
                                                  Slope_lower = Q2.5,
                                                  Slope_upper = Q97.5) %>% 
                                           dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(dataset_label) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'dataset_label')

Nstd_fragSize_group_coefs <- bind_cols(Nstd_fS_coef[[1]][,,'Intercept'] %>% 
                                          as_tibble() %>% 
                                          mutate(Intercept = Estimate,
                                                 Intercept_lower = Q2.5,
                                                 Intercept_upper = Q97.5,
                                                 dataset_label = rownames(Nstd_fS_coef[[1]][,,'Intercept'])) %>% 
                                          dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                       Nstd_fS_coef[[1]][,,'c.lfs'] %>% 
                                          as_tibble() %>% 
                                          mutate(Slope = Estimate,
                                                 Slope_lower = Q2.5,
                                                 Slope_upper = Q97.5) %>% 
                                          dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(dataset_label) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'dataset_label')

#---- regression plots showing study-level slopes-----
# setwd('~/Dropbox/Habitat loss meta-analysis/analysis/figs/')
setwd(paste(path2temp,"figs/", sep = ""))

# S_std_2 (Felix recommends this one for the main result)
S_std_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_std_2, colour = dataset_label),
             size = 1.5) +
  geom_segment(data = Sstd2_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = dataset_label,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Sstd2_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Sstd2_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.1, y = 132,
           label = paste("beta == ", #[Frag.~size]
                         round(Sstd2_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         " (",
                         round(Sstd2_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
                         " - ",
                         round(Sstd2_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
                         ")"),  
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[std]))) +
  theme_bw() +
  theme(legend.position = 'none')

# Sn
Sn_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_n, colour = dataset_label),
             size = 1.5) +
  geom_segment(data = Sn_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = dataset_label,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Sn_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Sn_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.1, y = 132,
           label = paste("beta == ", #[Frag.~size]
                         round(Sn_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         "  (",
                         round(Sn_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
                         " - ",
                         round(Sn_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
                         ")"),  
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[n]))) +
  theme_bw() +
  theme(legend.position = 'none')

# S_PIE
Spie_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_PIE, colour = dataset_label),
             size = 1.5) +
  geom_segment(data = S_PIE_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = dataset_label,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = S_PIE_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = S_PIE_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.1, y = 132,
           label = paste("beta == ", #[Frag.~size]
                         round(S_PIE_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         "  (",
                         round(S_PIE_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
                         " - ",
                         round(S_PIE_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
                         ")"),  
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[PIE]))) +
  theme_bw() +
  theme(legend.position = 'none')

# S_cov
Scov_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_cov, colour = dataset_label),
             size = 1.5) +
  geom_segment(data = Scov_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = dataset_label,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Scov_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Scov_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.1, y = 132,
           label = paste("beta == ", #[Frag.~size]
                         round(Scov_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         "  (",
                         round(Scov_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
                         " - ",
                         round(Scov_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
                         ")"),  
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[cov]))) +
  theme_bw() +
  theme(legend.position = 'none')

# S_cov
Schao_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = S_chao, colour = dataset_label),
             size = 1.5) +
  geom_segment(data = Schao_lognorm_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = dataset_label,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Schao_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Schao_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.1, y = 900,
           label = paste("beta == ", #[Frag.~size]
                         round(Schao_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         "  (",
                         round(Schao_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
                         " - ",
                         round(Schao_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
                         ")"),  
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(2, 16, 128, 1024)) +#
  scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[chao]))) +
  theme_bw() +
  theme(legend.position = 'none')

# N_std
Nstd_regPlot <- ggplot() +
  # data
  geom_point(data = frag,
             aes(x = frag_size_num, y = N_std, 
                 colour = dataset_label),
             size = 1.5) +
  geom_segment(data = Nstd_fragSize_group_coefs,
               aes(group = dataset_label,
                   colour = dataset_label,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Nstd_fS_fitted, 
            aes(x = frag_size_num,
                y = Estimate),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Nstd_fS_fitted,
              aes(x = frag_size_num,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.175, y = 20000,
           label = paste("beta == ", #[Frag.~size]
                         round(Nstd_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
                         "  (",
                         round(Nstd_lognorm_fragSize_fixef['c.lfs','Q2.5'],3),
                         " - ",
                         round(Nstd_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
                         ")"),  
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000)) +
  scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(N[std]))) +
  theme_bw() +
  theme(legend.position = 'none')

S_std_regPlot
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/S_std_fragSize.png', 
#        width = 200, height = 200, units = 'mm')
Sn_regPlot
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/S_n_fragSize.png', 
#        width = 200, height = 200, units = 'mm')
Scov_regPlot
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/S_cov_fragSize.png', 
#        width = 200, height = 200, units = 'mm')
Schao_regPlot
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/S_chao_fragSize.png', 
#        width = 200, height = 200, units = 'mm')
Spie_regPlot
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/S_PIE_fragSize.png', 
#        width = 200, height = 200, units = 'mm')

Nstd_regPlot
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/N_std_fragSize.png', 
#        width = 200, height = 200, units = 'mm')

# ggsave('~/fragmentation/figs/fragSize_regression.pdf', width = 250, height = 220, units = 'mm')
##---coef plots---------
# get the metadata...
meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep = ';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

# create unique dataframes with columns for plotting relationships between change in different metrics
# there are different studies retained for each metric, so put 'em together for each plot separately
S_std_study_slope <- Sstd2_lognorm_fragSize_group_coefs %>% 
    mutate(S_std2_slope = Slope,
           S_std2_lower = Slope_lower,
           S_std2_upper = Slope_upper) %>% 
    select(dataset_label, S_std2_slope, S_std2_lower, S_std2_upper)

Sn_study_slope<- Sn_lognorm_fragSize_group_coefs %>% 
    mutate(S_n_slope = Slope,
           S_n_lower = Slope_lower,
           S_n_upper = Slope_upper) %>% 
    select(dataset_label, S_n_slope, S_n_lower, S_n_upper)

Scov_study_slope<- Scov_lognorm_fragSize_group_coefs %>% 
  mutate(S_cov_slope = Slope,
         S_cov_lower = Slope_lower,
         S_cov_upper = Slope_upper) %>% 
  select(dataset_label, S_cov_slope, S_cov_lower, S_cov_upper)
