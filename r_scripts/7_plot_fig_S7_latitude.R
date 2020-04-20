## code to fit models of study-level richness slope as a function of absolute latitude, 
## with and without interaction between latitude and taxa
## models borrow from 'meta-analytic' modelling concept: they include uncertainty in the 
## study-level slope estimates

## finally, we'll plot results of best-fitting model

# get the coefficients for all the results
source(paste0(path2wd, 'r_scripts/5a_fragSize_coef_wrangle.R'))

meta <- read.csv(paste0(path2wd, 'data/new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id) %>% 
  separate(coordinates, into = c('y', 'x'), sep = ', ', remove = F) %>% 
  mutate(x = as.numeric(x),
         y = as.numeric(y),
         abs_lat = abs(y),
         latitude = climate,
         abs_lat_bin = cut(abs_lat,
                           breaks = seq(0,70, by = 10),
                           labels = c('0-10', '10-20','20-30','30-40',
                                      '40-50','50-60','60-70')))

S_std_study_slope <- Sstd2_lognorm_fragSize_group_coefs %>% 
  mutate(S_std_slope = Slope,
         error = error,
         S_std_lower = Slope_lower,
         S_std_upper = Slope_upper) %>% 
  select(dataset_label, S_std_slope, error, S_std_lower, S_std_upper) %>% 
  left_join(meta,
            by = 'dataset_label')


# include uncertainty of the study-level estimate in the model 
slope_lat <- brms::brm(bf(S_std_slope | se(error) ~ abs_lat + (1|dataset_label)),
                       data = S_std_study_slope,
                       cores = 4, chains = 4)

pp_check(slope_lat)
plot(slope_lat)

slope_lat_taxa <- brms::brm(bf(S_std_slope | se(error) ~ abs_lat*taxa + (1|dataset_label)),
                            data = S_std_study_slope,
                            cores = 4, chains = 4)

slope_lat <- add_criterion(slope_lat, 
                           criterion = c('loo', 'waic'))
slope_lat_taxa <- add_criterion(slope_lat_taxa, 
                                criterion = c('loo', 'waic'),
                                reloo = TRUE)

# not much support for including taxa (though the fits do not differ dramatically)
model_weights(slope_lat, slope_lat_taxa)

# model with taxa x latitude interaction wins: more than twice as likely
slope_lat_fixef <- cbind(slope_lat$data,
                         fitted(slope_lat, 
                                re_formula = NA)) %>% 
  as_tibble()

ggplot() +
  geom_point(data = S_std_study_slope,
             aes(x = abs(y), y = S_std_slope),
             size = 1) +
  geom_linerange(data = S_std_study_slope,
                 aes(x = abs(y), ymin = S_std_lower, ymax = S_std_upper),
                 size = 0.5, alpha = 0.5) +
  geom_line(data = slope_lat_fixef,
            aes(x = abs_lat, y = Estimate), 
            size = 1) +
  geom_ribbon(data = slope_lat_fixef,
              aes(x = abs_lat, ymin = Q2.5, ymax = Q97.5), 
              alpha = 0.5) +
  labs(x = '|Latitude|',
       y = 'Richness ~ fragment size slope estimate') +
  theme_bw() +
  theme(text = element_text(size = 7))

# 1.5 column figure size
# set local directory
# ggsave('~/Dropbox/Frag Database (new)/Manuscript for Nature/revision3/figures/Ex_Dat_Fig7.png',
#        width = 120, height = 120, units = 'mm')
