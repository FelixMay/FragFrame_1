# execute 0_init_dirs_load_packages.R for packages and directory

# load model fits
load(paste0(path2wd, 'main_results/Jtu_z1i_fragSize.Rdata'))
load(paste0(path2wd, 'main_results/Rtu_z1i_fragSize.Rdata'))

load(paste0(path2wd, 'main_results/Jne_zi_fragSize.Rdata'))
load(paste0(path2wd, 'main_results/Rne_zi_fragSize.Rdata'))

frag_beta <- read_csv(paste0(path2wd, 'intermediate_results/2_betapart_frag_fcont_10_mabund_as_is.csv'))

# get the metadata
meta <- read_delim(paste0(path2wd, 'data/new_meta_2_merge.csv'),  delim =';') %>% 
   dplyr::rename(dataset_label = dataset_id) 

frag_beta <- frag_beta %>% 
  # centre covariate before fitting
  mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))

frag_beta <- left_join(frag_beta, 
                       meta, by = 'dataset_label')

frag_beta$time.since.fragmentation <- factor(frag_beta$time.since.fragmentation,
                                                       levels = c('Recent (less than 20 years)',
                                                                  'Intermediate (20-100 years)',
                                                                  'long (100+ years)'),
                                                       labels = c('< 20 years',
                                                                  '20-100 years',
                                                                  '> 100 years'))

frag_beta$Matrix.category <- factor(frag_beta$Matrix.category,
                                              levels = c('light filter', 'intermediate', 'harsh filter'),
                                              labels = c('Light', 'Intermediate', 'Harsh'))

frag_beta$biome <- factor(frag_beta$biome,
                                    levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                                    labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))

frag_beta$taxa <- factor(frag_beta$taxa,
                                   levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                                   labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

frag_beta <- frag_beta %>% 
  unite(col = 'frag_matrix', c(sphere.fragment, sphere.matrix),
        remove = F, sep = ', ')


#------wrangle for plotting
# for plotting fixed effects----------------
Jtu_z1i_fitted <- cbind(Jtu_z1i_fragSize$data,
                         fitted(Jtu_z1i_fragSize, re_formula = NA, scale = 'linear')) %>%
  as_tibble() %>% 
  left_join(frag_beta %>% 
               distinct(dataset_label, cl10ra, frag_size_num.x, frag_size_num.y, log10_ratio_area),
             by = c('dataset_label', 'cl10ra'))# %>% 
  # want the 'fitted' values of the coi component too
  # bind_cols(
  #   fitted(Jtu_z1i_fragSize, re_formula = NA, dpar = 'coi', scale = 'linear') %>% 
  #     as_tibble() %>% 
  #     mutate(coi_Estimate = Estimate,
  #            coi_Q2.5 = Q2.5,
  #            coi_Q97.5 = Q97.5) %>% 
  #     select(coi_Estimate, coi_Q2.5, coi_Q97.5))


Jtu_z1i_fixef <- fixef(Jtu_z1i_fragSize)
# random effects
Jtu_z1i_coef <- coef(Jtu_z1i_fragSize)

Jtu_z1i_group_coefs <- bind_cols(Jtu_z1i_coef[[1]][,,'Intercept'] %>% 
                                                  as_tibble() %>% 
                                                  mutate(Intercept = Estimate,
                                                         Intercept_lower = Q2.5,
                                                         Intercept_upper = Q97.5,
                                                         dataset_label = rownames(Jtu_z1i_coef[[1]][,,'Intercept'])) %>% 
                                                  dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                 Jtu_z1i_coef[[1]][,,'cl10ra'] %>% 
                                                  as_tibble() %>% 
                                                  mutate(Slope = Estimate,
                                                         Slope_lower = Q2.5,
                                                         Slope_upper = Q97.5) %>% 
                                                  dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                 Jtu_z1i_coef[[1]][,,'zoi_cl10ra'] %>% 
                                   as_tibble() %>% 
                                   mutate(zoi_Slope = Estimate,
                                          zoi_Slope_lower = Q2.5,
                                          zoi_Slope_upper = Q97.5) %>% 
                                   dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag_beta %>% 
               group_by(dataset_label) %>% 
               summarise(xmin.x = min(frag_size_num.x),
                         xmax.x = max(frag_size_num.x),
                         xmin.y = min(frag_size_num.y),
                         xmax.y = max(frag_size_num.y),
                         cxmin = min(cl10ra),
                         cxmax = max(cl10ra)),
             by = 'dataset_label') %>% 
  left_join(meta, by = 'dataset_label')

Jtu_z1i_group_coefs$time.since.fragmentation <- factor(Jtu_z1i_group_coefs$time.since.fragmentation,
                                                 levels = c('Recent (less than 20 years)',
                                                            'Intermediate (20-100 years)',
                                                            'long (100+ years)'),
                                                 labels = c('< 20 years',
                                                            '20-100 years',
                                                            '> 100 years'))

Jtu_z1i_group_coefs$Matrix.category <- factor(Jtu_z1i_group_coefs$Matrix.category,
                                        levels = c('light filter', 'intermediate', 'harsh filter'),
                                        labels = c('Light', 'Intermediate', 'Harsh'))

Jtu_z1i_group_coefs$biome <- factor(Jtu_z1i_group_coefs$biome,
                              levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                              labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))

Jtu_z1i_group_coefs$taxa <- factor(Jtu_z1i_group_coefs$taxa,
                             levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                             labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

Jtu_z1i_group_coefs <- Jtu_z1i_group_coefs %>% 
  unite(col = 'frag_matrix', c(sphere.fragment, sphere.matrix),
        remove = F, sep = ', ')


# repeat for Ruzicka
Rtu_z1i_fitted <- cbind(Rtu_z1i_fragSize$data,
                        fitted(Rtu_z1i_fragSize, re_formula = NA, scale = 'linear')) %>%
  as_tibble() %>% 
  left_join(frag_beta %>% 
              distinct(dataset_label, cl10ra, frag_size_num.x, frag_size_num.y, log10_ratio_area),
            by = c('dataset_label', 'cl10ra')) #%>% 
  # want the 'fitted' values of the coi component too
  # bind_cols(
  #   fitted(Rtu_z1i_fragSize, re_formula = NA, dpar = 'coi', scale = 'linear') %>% 
  #     as_tibble() %>% 
  #     mutate(coi_Estimate = Estimate,
  #            coi_Q2.5 = Q2.5,
  #            coi_Q97.5 = Q97.5) %>% 
  #     select(coi_Estimate, coi_Q2.5, coi_Q97.5))


Rtu_z1i_fixef <- fixef(Rtu_z1i_fragSize)

Rtu_z1i_coef <- coef(Rtu_z1i_fragSize)

Rtu_z1i_group_coefs <- bind_cols(Rtu_z1i_coef[[1]][,,'Intercept'] %>% 
                                   as_tibble() %>% 
                                   mutate(Intercept = Estimate,
                                          Intercept_lower = Q2.5,
                                          Intercept_upper = Q97.5,
                                          dataset_label = rownames(Rtu_z1i_coef[[1]][,,'Intercept'])) %>% 
                                   dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                 Rtu_z1i_coef[[1]][,,'cl10ra'] %>% 
                                   as_tibble() %>% 
                                   mutate(Slope = Estimate,
                                          Slope_lower = Q2.5,
                                          Slope_upper = Q97.5) %>% 
                                   dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag_beta %>% 
               group_by(dataset_label) %>% 
               summarise(xmin.x = min(frag_size_num.x),
                         xmax.x = max(frag_size_num.x),
                         xmin.y = min(frag_size_num.y),
                         xmax.y = max(frag_size_num.y),
                         cxmin = min(cl10ra),
                         cxmax = max(cl10ra)),
             by = 'dataset_label') %>% 
  left_join(meta, by = 'dataset_label')

Rtu_z1i_group_coefs$time.since.fragmentation <- factor(Rtu_z1i_group_coefs$time.since.fragmentation,
                                                       levels = c('Recent (less than 20 years)',
                                                                  'Intermediate (20-100 years)',
                                                                  'long (100+ years)'),
                                                       labels = c('< 20 years',
                                                                  '20-100 years',
                                                                  '> 100 years'))

Rtu_z1i_group_coefs$Matrix.category <- factor(Rtu_z1i_group_coefs$Matrix.category,
                                              levels = c('light filter', 'intermediate', 'harsh filter'),
                                              labels = c('Light', 'Intermediate', 'Harsh'))

Rtu_z1i_group_coefs$biome <- factor(Rtu_z1i_group_coefs$biome,
                                    levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                                    labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))

Rtu_z1i_group_coefs$taxa <- factor(Rtu_z1i_group_coefs$taxa,
                                   levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                                   labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

Rtu_z1i_group_coefs <- Rtu_z1i_group_coefs %>% 
  unite(col = 'frag_matrix', c(sphere.fragment, sphere.matrix),
        remove = F, sep = ', ')


# now for the nestedness components
Jne_zi_fitted <- cbind(Jne_zi_fragSize$data,
                        fitted(Jne_zi_fragSize, re_formula = NA, scale = 'linear')) %>%
  as_tibble() %>% 
  left_join(frag_beta %>%
              distinct(dataset_label, cl10ra, frag_size_num.x, frag_size_num.y, log10_ratio_area),
            by = c('dataset_label', 'cl10ra')) #%>% 
  # want the 'fitted' values of the coi component too
  # bind_cols(
  #   fitted(Jne_zi_fragSize, re_formula = NA, dpar = 'zi', scale = 'linear') %>% 
  #     as_tibble() %>% 
  #     mutate(zi_Estimate = Estimate,
  #            zi_Q2.5 = Q2.5,
  #            zi_Q97.5 = Q97.5) %>% 
  #     select(zi_Estimate, zi_Q2.5, zi_Q97.5))


Jne_zi_fixef <- fixef(Jne_zi_fragSize)
# random effects
Jne_zi_coef <- coef(Jne_zi_fragSize)

Jne_zi_group_coefs <- bind_cols(Jne_zi_coef[[1]][,,'Intercept'] %>% 
                                   as_tibble() %>% 
                                   mutate(Intercept = Estimate,
                                          Intercept_lower = Q2.5,
                                          Intercept_upper = Q97.5,
                                          dataset_label = rownames(Jne_zi_coef[[1]][,,'Intercept'])) %>% 
                                   dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                 Jne_zi_coef[[1]][,,'cl10ra'] %>% 
                                   as_tibble() %>% 
                                   mutate(Slope = Estimate,
                                          Slope_lower = Q2.5,
                                          Slope_upper = Q97.5) %>% 
                                   dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag_beta %>% 
               filter(method == 'Baselga family, Jaccard') %>% 
               group_by(dataset_label) %>% 
               summarise(xmin.x = min(frag_size_num.x),
                         xmax.x = max(frag_size_num.x),
                         xmin.y = min(frag_size_num.y),
                         xmax.y = max(frag_size_num.y),
                         cxmin = min(cl10ra),
                         cxmax = max(cl10ra)),
             by = 'dataset_label') %>% 
  left_join(meta, by = 'dataset_label')

Jne_zi_group_coefs$time.since.fragmentation <- factor(Jne_zi_group_coefs$time.since.fragmentation,
                                                       levels = c('Recent (less than 20 years)',
                                                                  'Intermediate (20-100 years)',
                                                                  'long (100+ years)'),
                                                       labels = c('< 20 years',
                                                                  '20-100 years',
                                                                  '> 100 years'))

Jne_zi_group_coefs$Matrix.category <- factor(Jne_zi_group_coefs$Matrix.category,
                                              levels = c('light filter', 'intermediate', 'harsh filter'),
                                              labels = c('Light', 'Intermediate', 'Harsh'))

Jne_zi_group_coefs$biome <- factor(Jne_zi_group_coefs$biome,
                                    levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                                    labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))

Jne_zi_group_coefs$taxa <- factor(Jne_zi_group_coefs$taxa,
                                   levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                                   labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

Jne_zi_group_coefs <- Jne_zi_group_coefs %>% 
  unite(col = 'frag_matrix', c(sphere.fragment, sphere.matrix),
        remove = F, sep = ', ')


# repeat for Ruzicka
Rne_zi_fitted <- cbind(Rne_zi_fragSize$data,
                        fitted(Rne_zi_fragSize, re_formula = NA, scale = 'linear')) %>%
  as_tibble() %>% 
  left_join(frag_beta %>% 
              distinct(dataset_label, cl10ra, frag_size_num.x, frag_size_num.y, log10_ratio_area),
            by = c('dataset_label', 'cl10ra')) #%>% 
  # want the 'fitted' values of the coi component too
  # bind_cols(
  #   fitted(Rne_zi_fragSize, re_formula = NA, dpar = 'zi', scale = 'linear') %>% 
  #     as_tibble() %>% 
  #     mutate(zi_Estimate = Estimate,
  #            zi_Q2.5 = Q2.5,
  #            zi_Q97.5 = Q97.5) %>% 
  #     select(zi_Estimate, zi_Q2.5, zi_Q97.5))


Rne_zi_fixef <- fixef(Rne_zi_fragSize)

Rne_zi_coef <- coef(Rne_zi_fragSize)

Rne_zi_group_coefs <- bind_cols(Rne_zi_coef[[1]][,,'Intercept'] %>% 
                                   as_tibble() %>% 
                                   mutate(Intercept = Estimate,
                                          Intercept_lower = Q2.5,
                                          Intercept_upper = Q97.5,
                                          dataset_label = rownames(Rne_zi_coef[[1]][,,'Intercept'])) %>% 
                                   dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                 Rne_zi_coef[[1]][,,'cl10ra'] %>% 
                                   as_tibble() %>% 
                                   mutate(Slope = Estimate,
                                          Slope_lower = Q2.5,
                                          Slope_upper = Q97.5) %>% 
                                   dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag_beta %>% 
               group_by(dataset_label) %>% 
               summarise(xmin.x = min(frag_size_num.x),
                         xmax.x = max(frag_size_num.x),
                         xmin.y = min(frag_size_num.y),
                         xmax.y = max(frag_size_num.y),
                         cxmin = min(cl10ra),
                         cxmax = max(cl10ra)),
             by = 'dataset_label') %>% 
  left_join(meta, by = 'dataset_label')

Rne_zi_group_coefs$time.since.fragmentation <- factor(Rne_zi_group_coefs$time.since.fragmentation,
                                                       levels = c('Recent (less than 20 years)',
                                                                  'Intermediate (20-100 years)',
                                                                  'long (100+ years)'),
                                                       labels = c('< 20 years',
                                                                  '20-100 years',
                                                                  '> 100 years'))

Rne_zi_group_coefs$Matrix.category <- factor(Rne_zi_group_coefs$Matrix.category,
                                              levels = c('light filter', 'intermediate', 'harsh filter'),
                                              labels = c('Light', 'Intermediate', 'Harsh'))

Rne_zi_group_coefs$biome <- factor(Rne_zi_group_coefs$biome,
                                    levels = c('forest', 'grassland', 'shrubland/steppe', 'wetland'),
                                    labels = c('Forest', 'Grassland', 'Shrubland or steppe', 'Wetland'))

Rne_zi_group_coefs$taxa <- factor(Rne_zi_group_coefs$taxa,
                                   levels = c('amphibians & reptiles', 'birds', 'invertebrates', 'mammals', 'plants'),
                                   labels = c('Amphibians & reptiles', 'Birds', 'Invertebrates', 'Mammals', 'Plants'))

Rne_zi_group_coefs <- Rne_zi_group_coefs %>% 
  unite(col = 'frag_matrix', c(sphere.fragment, sphere.matrix),
        remove = F, sep = ', ')


