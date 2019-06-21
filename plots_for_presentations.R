# code to plot regressions for all metrics (no interaction models)
library(tidyverse)
library(brms)
# source('~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/fragSize_coef_wrangle.R')

load(path2wd %+% "fragSize_all_coefs.RData")

#---- regression plots showing study-level slopes-----
setwd(paste(path2Dropbox,"analysis_apr19/figs_presentation_felix/", sep = ""))


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
  # # add regression coefficient and uncertainty interval
  # annotate('text', x = 0.1, y = 132,
  #          label = paste("beta == ", #[Frag.~size]
  #                        round(Sstd2_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
  #                        " (",
  #                        round(Sstd2_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
  #                        " - ",
  #                        round(Sstd2_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
  #                        ")"),
  #          parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[std]))) +
  theme_bw() +
  theme(legend.position = 'none')
# 
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
  # # add regression coefficient and uncertainty interval
  # annotate('text', x = 0.1, y = 132,
  #          label = paste("beta == ", #[Frag.~size]
  #                        round(Sn_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
  #                        "  (",
  #                        round(Sn_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
  #                        " - ",
  #                        round(Sn_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
  #                        ")"),
  #          parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
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
  # # add regression coefficient and uncertainty interval
  # annotate('text', x = 0.1, y = 132,
  #          label = paste("beta == ", #[Frag.~size]
  #                        round(S_PIE_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
  #                        "  (",
  #                        round(S_PIE_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
  #                        " - ",
  #                        round(S_PIE_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
  #                        ")"),
  #          parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[PIE]))) +
  theme_bw() +
  theme(legend.position = 'none')

# # S_cov
# Scov_regPlot <- ggplot() +
#   # data
#   geom_point(data = frag,
#              aes(x = frag_size_num, y = S_cov, colour = dataset_label),
#              size = 1.5) +
#   geom_segment(data = Scov_lognorm_fragSize_group_coefs,
#                aes(group = dataset_label,
#                    colour = dataset_label,
#                    x = xmin,
#                    xend = xmax,
#                    y = exp(Intercept + Slope * cxmin),
#                    yend = exp(Intercept + Slope * cxmax)),
#                size = 0.5) +
#   # fixed effect
#   geom_line(data = Scov_fS_fitted, 
#             aes(x = frag_size_num,
#                 y = Estimate),
#             size = 1.5) +
#   # fixed effect uncertainty
#   geom_ribbon(data = Scov_fS_fitted,
#               aes(x = frag_size_num,
#                   ymin = Q2.5,
#                   ymax = Q97.5),
#               alpha = 0.3) +
#   # add regression coefficient and uncertainty interval
#   annotate('text', x = 0.1, y = 132,
#            label = paste("beta == ", #[Frag.~size]
#                          round(Scov_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
#                          "  (",
#                          round(Scov_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
#                          " - ",
#                          round(Scov_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
#                          ")"),  
#            parse = T) +
#   scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#   scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
#   scale_colour_viridis_d(guide=F) +
#   labs(x = 'Fragment size (hectares)',
#        y = expression(paste(S[cov]))) +
#   theme_bw() +
#   theme(legend.position = 'none')
# 
# # S_cov
# Schao_regPlot <- ggplot() +
#   # data
#   geom_point(data = frag,
#              aes(x = frag_size_num, y = S_chao, colour = dataset_label),
#              size = 1.5) +
#   geom_segment(data = Schao_lognorm_fragSize_group_coefs,
#                aes(group = dataset_label,
#                    colour = dataset_label,
#                    x = xmin,
#                    xend = xmax,
#                    y = exp(Intercept + Slope * cxmin),
#                    yend = exp(Intercept + Slope * cxmax)),
#                size = 0.5) +
#   # fixed effect
#   geom_line(data = Schao_fS_fitted, 
#             aes(x = frag_size_num,
#                 y = Estimate),
#             size = 1.5) +
#   # fixed effect uncertainty
#   geom_ribbon(data = Schao_fS_fitted,
#               aes(x = frag_size_num,
#                   ymin = Q2.5,
#                   ymax = Q97.5),
#               alpha = 0.3) +
#   # add regression coefficient and uncertainty interval
#   annotate('text', x = 0.1, y = 900,
#            label = paste("beta == ", #[Frag.~size]
#                          round(Schao_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
#                          "  (",
#                          round(Schao_lognorm_fragSize_fixef['c.lfs','Q2.5'],2),
#                          " - ",
#                          round(Schao_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
#                          ")"),  
#            parse = T) +
#   scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#   scale_y_continuous(trans = 'log', breaks = c(2, 16, 128, 1024)) +#
#   scale_colour_viridis_d(guide=F) +
#   labs(x = 'Fragment size (hectares)',
#        y = expression(paste(S[chao]))) +
#   theme_bw() +
#   theme(legend.position = 'none')

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
  # annotate('text', x = 0.01, y = 5000, hjust = 0, vjust = 0,
  #          label = paste("beta == ", #[Frag.~size]
  #                        round(Nstd_lognorm_fragSize_fixef['c.lfs','Estimate'],2),
  #                        "  (",
  #                        round(Nstd_lognorm_fragSize_fixef['c.lfs','Q2.5'],3),
  #                        " - ",
  #                        round(Nstd_lognorm_fragSize_fixef['c.lfs','Q97.5'],2),
  #                        ")"),  
  #          parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000)) +
  # scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(N[std]))) +
  theme_bw() +
  theme(legend.position = 'none')

# S_std_regPlot
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/S_std_fragSize.png',
#        width = 200, height = 200, units = 'mm')

# Sn_regPlot
ggsave("S_n_fragSize.png", Sn_regPlot, width = 110, height = 100, units = 'mm')

# Scov_regPlot
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/S_cov_fragSize.png',
#        width = 200, height = 200, units = 'mm')
# Schao_regPlot
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/S_chao_fragSize.png',
#        width = 200, height = 200, units = 'mm')

# Spie_regPlot
ggsave("S_PIE_fragSize.png", Spie_regPlot, width = 110, height = 100, units = 'mm')

# Nstd_regPlot
ggsave("N_std_fragSize.png", Nstd_regPlot, width = 110, height = 100, units = 'mm')

