# code to plot beta-diversity regression
rm(list=ls())
library(tidyverse)
library(brms)
source('~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/beta_frag_coef_wrangle.R')

##--------get study-level sd for adding to figure----
study_sd <- bind_rows(
  cbind.data.frame(
    model = 'Jtu',
    Jtu_z1i_fS %>% 
  tidybayes::spread_draws(sd_dataset_label__cl10ra) %>% 
  tidybayes::mean_qi()),
  cbind.data.frame(
    model = 'Rtu',
    Rtu_z1i_fS %>% 
      tidybayes::spread_draws(sd_dataset_label__cl10ra) %>% 
      tidybayes::mean_qi()),
  cbind.data.frame(
    model = 'Jne',
    Jne_zi_fragSize %>% 
      tidybayes::spread_draws(sd_dataset_label__cl10ra) %>% 
      tidybayes::mean_qi()),
  cbind.data.frame(
    model = 'Rne',
    Rne_zi_fragSize %>% 
      tidybayes::spread_draws(sd_dataset_label__cl10ra) %>% 
      tidybayes::mean_qi()))
  
##--------now plot the regressions----------------
turnover_regression <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Jaccard'),
             aes(x = frag_size_num.x/frag_size_num.y, y = repl, colour = taxa),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Jtu_z1i_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
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
  # fixed effect uncertainty
  geom_ribbon(data = Jtu_z1i_fitted,
              aes(x = frag_size_num.x/frag_size_num.y,
                  ymin = plogis(Q2.5),
                  ymax = plogis(Q97.5)),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  # annotate('text', x = Inf, y = 0.1, hjust = 1.2, vjust = 0,
  #        label = paste("beta[Jtu] == ", #[Frag.~size]
  #                      round(Jtu_z1i_fixef['cl10ra','Estimate'],2),
  #                      "  (",
  #                      round(Jtu_z1i_fixef['cl10ra','Q2.5'],2),
  #                      " - ",
  #                      round(Jtu_z1i_fixef['cl10ra','Q97.5'],2),
  #                      ")"),
  #        parse = T) +
  # annotate('text', x = Inf, y = 0, hjust = 1.2, vjust = 0,
  #          label = paste("sigma[beta[Jtu]*','~study] == ", #[Frag.~size]
  #                        round(study_sd %>% filter(model=='Jtu') %>% 
  #                                .$sd_dataset_label__cl10ra,2),
  #                        "  (",
  #                        round(study_sd %>% filter(model=='Jtu') %>% 
  #                                .$.lower,2),
  #                        " - ",
  #                        round(study_sd %>% filter(model=='Jtu') %>% 
  #                                .$.upper,2),
  #                        ")"),
  #          parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = 'Turnover (Jaccard)',
       tag = 'a') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 13))


turnover_regression_ruzicka <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Ruzicka'),
             aes(x = frag_size_num.x/frag_size_num.y, y = repl, colour = taxa),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Rtu_z1i_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
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
  # annotate('text', x = Inf, y = 0.1, hjust = 1.2, vjust = 0,
  #          label = paste("beta[Rbal] == ", #[Frag.~size]
  #                        round(Rtu_z1i_fixef['cl10ra','Estimate'],2),
  #                        "  (",
  #                        round(Rtu_z1i_fixef['cl10ra','Q2.5'],2),
  #                        " - ",
  #                        round(Rtu_z1i_fixef['cl10ra','Q97.5'],2),
  #                        ")"),
  #          parse = T) +
  # annotate('text', x = Inf, y = 0, hjust = 1.2, vjust = 0,
  #          label = paste("sigma[study, fS] == ", #[Frag.~size]
  #                        round(study_sd %>% filter(model=='Rtu') %>% 
  #                                .$sd_dataset_label__cl10ra,2),
  #                        "  (",
  #                        round(study_sd %>% filter(model=='Rtu') %>% 
  #                                .$.lower,2),
  #                        " - ",
  #                        round(study_sd %>% filter(model=='Rtu') %>% 
  #                                .$.upper,2),
  #                        ")"),
  #          parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = 'Balanced abundance (Ruzicka)',
       tag = 'b') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 13))

nestedness_regression <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Jaccard'),
             aes(x = frag_size_num.x/frag_size_num.y, y = rich, colour = taxa),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Jne_zi_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
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
  # annotate('text', x = Inf, y = 1, hjust = 1.2, vjust = 1.4,
  #        label = paste("beta[Jne] == ", #[Frag.~size]
  #                      round(Jne_zi_fixef['cl10ra','Estimate'],2),
  #                      "  (",
  #                      round(Jne_zi_fixef['cl10ra','Q2.5'],2),
  #                      " - ",
  #                      round(Jne_zi_fixef['cl10ra','Q97.5'],2),
  #                      ")"),
  #        parse = T) +
  # annotate('text', x = Inf, y = 0.8, hjust = 1.2, vjust = 0,
  #          label = paste("sigma[study, fS] == ", #[Frag.~size]
  #                        round(study_sd %>% filter(model=='Jne') %>% 
  #                                .$sd_dataset_label__cl10ra,2),
  #                        "  (",
  #                        round(study_sd %>% filter(model=='Jne') %>% 
  #                                .$.lower,2),
  #                        " - ",
  #                        round(study_sd %>% filter(model=='Jne') %>% 
  #                                .$.upper,2),
  #                        ")"),
  #          parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = 'Nestedness (Jaccard)', 
       tag = 'c') +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(legend.position = 'none', text = element_text(size = 13))

nestedness_regression_ruzicka <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Ruzicka'),
             aes(x = frag_size_num.x/frag_size_num.y, y = rich, colour = taxa),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Rne_zi_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
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
  # annotate('text', x = Inf, y = 1, hjust = 1.2, vjust = 1.4,
  #          label = paste("beta[Rgrad] == ", #[Frag.~size]
  #                        round(Rne_zi_fixef['cl10ra','Estimate'],2),
  #                        "  (",
  #                        round(Rne_zi_fixef['cl10ra','Q2.5'],2),
  #                        " - ",
  #                        round(Rne_zi_fixef['cl10ra','Q97.5'],2),
  #                        ")"),
  #          parse = T) +
  # annotate('text', x = Inf, y = 0.8, hjust = 1.2, vjust = 0,
  #          label = paste("sigma[study, fS] == ", #[Frag.~size]
  #                        round(study_sd %>% filter(model=='Rne') %>% 
  #                                .$sd_dataset_label__cl10ra,2),
  #                        "  (",
  #                        round(study_sd %>% filter(model=='Rne') %>% 
  #                                .$.lower,2),
  #                        " - ",
  #                        round(study_sd %>% filter(model=='Rne') %>% 
  #                                .$.upper,2),
  #                        ")"),
  #          parse = T) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  labs(x = '',
       y = 'Abundance gradient (Ruzicka)',
       tag = 'd') +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 13))

# plot to generate legend 
taxa_legend <- ggplot() +
  # data
  geom_point(data = frag_beta %>% filter(method == 'Baselga family, Ruzicka'),
             aes(x = frag_size_num.x/frag_size_num.y, y = rich, colour = taxa),
             size = 0.5,
             alpha = 0.15) +
  geom_segment(data = Rne_zi_group_coefs,
               aes(group = dataset_label,
                   colour = taxa,
                   x = xmin.x/xmin.y,
                   xend = xmax.x/xmin.y,
                   y = plogis(Intercept + Slope * cxmin),
                   yend = plogis(Intercept + Slope * cxmax)),
               size = 0.5) +
  scale_colour_brewer(name = 'Taxa', type = 'qual', palette = 'Dark2') +
  theme_bw() +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        text = element_text(size = 13)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))

source('~/Dropbox/1current/R_random/functions/gg_legend.R')
taxa_colour = gg_legend(taxa_legend)

top <- cowplot::plot_grid(taxa_colour)
bottom <- cowplot::plot_grid(turnover_regression,
                             turnover_regression_ruzicka,
                             nestedness_regression,
                             nestedness_regression_ruzicka,
                   nrow = 2, align = 'hv')

cowplot::plot_grid(top, bottom, 
                   nrow = 2,
                   rel_heights = c(0.075,1)) +
  cowplot::draw_label('Ratio of fragment sizes (log-scale)',
                      y = 0.02)

# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/fig3_grey.png', width = 250, height = 220, units = 'mm')
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/figS4_taxa_colour.png', width = 250, height = 220, units = 'mm')
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/beta_frag.png', width = 250, height = 125, units = 'mm')