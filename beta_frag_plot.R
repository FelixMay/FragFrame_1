# code to plot beta-diversity regression
rm(list=ls())
library(tidyverse)
library(brms)
source('~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/beta_frag_coef_wrangle.R')

##--------now plot the regressions----------------
turnover_regression <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Jaccard'),
             aes(x = frag_size_num.x/frag_size_num.y, y = repl, colour = dataset_label),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Jtu_z1i_group_coefs,
               aes(group = dataset_label,
                   colour = dataset_label,
                   x = xmin.x/xmin.y,
                   # trick: max x-value requires a min in the denominator 
                   xend = xmax.x/xmin.y,
                   y = plogis(Intercept + Slope * cxmin),
                   yend = plogis(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Jtu_z1i_fitted, 
            aes(x = frag_size_num.x/frag_size_num.y,
                y = plogis(Estimate)),
            size = 1.5) +
  # fixed effect for Ruzicka
  # geom_line(data = Rtu_z1i_fitted,
  #           aes(x = frag_size_num.x/frag_size_num.y,
  #               y = plogis(Estimate)),
  #           size = 1, lty = 3) +
  # fixed effect uncertainty
  geom_ribbon(data = Jtu_z1i_fitted,
              aes(x = frag_size_num.x/frag_size_num.y,
                  ymin = plogis(Q2.5),
                  ymax = plogis(Q97.5)),
              alpha = 0.3) +
  # fixed effect uncertainty: Ruzicka
  # geom_ribbon(data = Rtu_z1i_fitted,
  #             aes(x = frag_size_num.x/frag_size_num.y,
  #                 ymin = plogis(Q2.5),
  #                 ymax = plogis(Q97.5)),
  #             alpha = 0.3) +
  # conditional probability of being one (given bernoulli mix)
  # geom_line(data = Jtu_z1i_fitted, 
  #           aes(x = frag_size_num.x/frag_size_num.y,
  #               y = plogis(coi_Estimate)),
  #           size = 1.5) +
  # geom_ribbon(data = Jtu_z1i_fitted,
  #             aes(x = frag_size_num.x/frag_size_num.y,
  #                 ymin = plogis(coi_Q2.5),
  #                 ymax = plogis(coi_Q97.5)),
  #             alpha = 0.3) +
  # add regression coefficient and uncertainty interval
annotate('text', x = Inf, y = 0, hjust = 1.2, vjust = 0,
         label = paste("beta[Jtu] == ", #[Frag.~size]
                       round(Jtu_z1i_fixef['cl10ra','Estimate'],2),
                       "  (",
                       round(Jtu_z1i_fixef['cl10ra','Q2.5'],2),
                       " - ",
                       round(Jtu_z1i_fixef['cl10ra','Q97.5'],2),
                       ")"),
         parse = T) +
  # annotate('text', x = 10^3.5, y = 0.9,
  #          label = paste("beta[Rbal] == ", #[Frag.~size]
  #                        round(Rtu_z1i_fixef['cl10ra','Estimate'],2),
  #                        "  (",
  #                        round(Rtu_z1i_fixef['cl10ra','Q2.5'],2),
  #                        " - ",
  #                        round(Rtu_z1i_fixef['cl10ra','Q97.5'],2),
  #                        ")"),
  #          parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  scale_color_grey(guide = F) +
  labs(x = 'Ratio of fragment sizes (log-scale)',
       y = 'Turnover (Jaccard)',
       tag = 'a') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 13))

turnover_regression_ruzicka <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Ruzicka'),
             aes(x = frag_size_num.x/frag_size_num.y, y = repl, colour = dataset_label),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Rtu_z1i_group_coefs,
               aes(group = dataset_label,
                   colour = dataset_label,
                   x = xmin.x/xmin.y,
                   # trick: max x-value requires a min in the denominator 
                   xend = xmax.x/xmin.y,
                   y = plogis(Intercept + Slope * cxmin),
                   yend = plogis(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect for Ruzicka
  geom_line(data = Rtu_z1i_fitted,
            aes(x = frag_size_num.x/frag_size_num.y,
                y = plogis(Estimate)),
            size = 1, lty = 1) +
  # fixed effect uncertainty: Ruzicka
  geom_ribbon(data = Rtu_z1i_fitted,
              aes(x = frag_size_num.x/frag_size_num.y,
                  ymin = plogis(Q2.5),
                  ymax = plogis(Q97.5)),
              alpha = 0.3) +
  annotate('text', x = Inf, y = 0, hjust = 1.2, vjust = 0,
           label = paste("beta[Rbal] == ", #[Frag.~size]
                         round(Rtu_z1i_fixef['cl10ra','Estimate'],2),
                         "  (",
                         round(Rtu_z1i_fixef['cl10ra','Q2.5'],2),
                         " - ",
                         round(Rtu_z1i_fixef['cl10ra','Q97.5'],2),
                         ")"),
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  scale_color_grey(guide = F) +
  labs(x = 'Ratio of fragment sizes (log-scale)',
       y = 'Balanced abundance (Ruzicka)',
       tag = 'c') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 13))

nestedness_regression <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Jaccard'),
             aes(x = frag_size_num.x/frag_size_num.y, y = rich, colour = dataset_label),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Jne_zi_group_coefs,
               aes(group = dataset_label,
                   colour = dataset_label,
                   x = xmin.x/xmin.y,
                   xend = xmax.x/xmin.y,
                   y = plogis(Intercept + Slope * cxmin),
                   yend = plogis(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Jne_zi_fitted, 
            aes(x = frag_size_num.x/frag_size_num.y,
                y = plogis(Estimate)),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = Jne_zi_fitted,
              aes(x = frag_size_num.x/frag_size_num.y,
                  ymin = plogis(Q2.5),
                  ymax = plogis(Q97.5)),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = Inf, y = 1, hjust = 1.2, vjust = 1.4,
         label = paste("beta[Jne] == ", #[Frag.~size]
                       round(Jne_zi_fixef['cl10ra','Estimate'],2),
                       "  (",
                       round(Jne_zi_fixef['cl10ra','Q2.5'],2),
                       " - ",
                       round(Jne_zi_fixef['cl10ra','Q97.5'],2),
                       ")"),
         parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  scale_colour_grey(guide = F) +
  labs(x = 'Ratio of fragment sizes (log-scale)',
       y = 'Nestedness (Jaccard)') +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(legend.position = 'none', text = element_text(size = 13))

nestedness_regression_ruzicka <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Ruzicka'),
             aes(x = frag_size_num.x/frag_size_num.y, y = rich, colour = dataset_label),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Rne_zi_group_coefs,
               aes(group = dataset_label,
                   colour = dataset_label,
                   x = xmin.x/xmin.y,
                   xend = xmax.x/xmin.y,
                   y = plogis(Intercept + Slope * cxmin),
                   yend = plogis(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Rne_zi_fitted,
            aes(x = frag_size_num.x/frag_size_num.y,
                y = plogis(Estimate)),
            size = 1, linetype = 1) +
  # fixed effect uncertainty: Ruzicka
  geom_ribbon(data = Rne_zi_fitted,
              aes(x = frag_size_num.x/frag_size_num.y,
                  ymin = plogis(Q2.5),
                  ymax = plogis(Q97.5)),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = Inf, y = 1, hjust = 1.2, vjust = 1.4,
           label = paste("beta[Rgrad] == ", #[Frag.~size]
                         round(Rne_zi_fixef['cl10ra','Estimate'],2),
                         "  (",
                         round(Rne_zi_fixef['cl10ra','Q2.5'],2),
                         " - ",
                         round(Rne_zi_fixef['cl10ra','Q97.5'],2),
                         ")"),
           parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  scale_colour_grey(guide = F) +
  labs(x = 'Ratio of fragment sizes (log-scale)',
       y = 'Abundance gradient (Ruzicka)',
       tag = 'd') +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 13))

cowplot::plot_grid(turnover_regression,
                   nestedness_regression,
                   turnover_regression_ruzicka,
                   nestedness_regression_ruzicka,
                   nrow = 2, align = 'hv')

ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/fig3_grey.png', width = 250, height = 220, units = 'mm')
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/beta_frag.png', width = 250, height = 125, units = 'mm')