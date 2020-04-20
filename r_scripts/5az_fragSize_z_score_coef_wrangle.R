# need to execute 0_init_dirs_load_packages.R first

# code to wrangle the coefficients for z-score regressions
load(paste0(path2wd, 'intermediate_results/fragSize_z_score_ref.Rdata'))

frag <- read_csv(paste0(path2wd, 'intermediate_results/2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

# load the meta data
meta <- read.csv(paste0(path2wd, 'data/new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')

#------wrangle for plotting
# for plotting fixed effects----------------
z_Sstd_fS_fitted <- cbind(z_Sstd_studT_fragSize$data,
                        fitted(z_Sstd_studT_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, frag_size_num),
             by = c('dataset_label', 'c.lfs'))

z_Sstd_studT_fragSize_fixef <- fixef(z_Sstd_studT_fragSize)

z_Sn_fS_fitted <- cbind(z_Sn_studT_fragSize$data,
                      fitted(z_Sn_studT_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, frag_size_num),
             by = c('dataset_label', 'c.lfs'))

z_Sn_studT_fragSize_fixef <- fixef(z_Sn_studT_fragSize)

z_Scov_fS_fitted <- cbind(z_Scov_studT_fragSize$data,
                        fitted(z_Scov_studT_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, frag_size_num),
             by = c('dataset_label', 'c.lfs'))

z_Scov_studT_fragSize_fixef <- fixef(z_Scov_studT_fragSize)

z_Schao_fS_fitted <- cbind(z_S_chao_studT_fragSize$data,
                         fitted(z_S_chao_studT_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, frag_size_num),
             by = c('dataset_label', 'c.lfs'))

z_Schao_lognorm_fragSize_fixef <- fixef(z_S_chao_studT_fragSize)

z_S_PIE_fS_fitted <- cbind(z_S_PIE_studT_fragSize$data,
                         fitted(z_S_PIE_studT_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(dataset_label, c.lfs, frag_size_num),
             by = c('dataset_label', 'c.lfs'))

z_S_PIE_studT_fragSize_fixef <- fixef(z_S_PIE_studT_fragSize)


# for plotting the random-effects----------------
z_Sstd_studT_fragSize_coef <- coef(z_Sstd_studT_fragSize)
z_Sn_studT_fragSize_coef <- coef(z_Sn_studT_fragSize)
z_Scov_studT_fragSize_coef <- coef(z_Scov_studT_fragSize)
z_Schao_studT_fragSize_coef <- coef(z_S_chao_studT_fragSize)
z_S_PIE_fS_coef <- coef(z_S_PIE_studT_fragSize)


z_Sstd_studT_fragSize_group_coefs <- bind_cols(z_Sstd_studT_fragSize_coef[[1]][,,'Intercept'] %>% 
                                                 as_tibble() %>% 
                                                 mutate(Intercept = Estimate,
                                                        Intercept_lower = Q2.5,
                                                        Intercept_upper = Q97.5,
                                                        dataset_label = rownames(z_Sstd_studT_fragSize_coef[[1]][,,'Intercept'])) %>% 
                                                 dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                               z_Sstd_studT_fragSize_coef[[1]][,,'c.lfs'] %>% 
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
                         cxmax = max(c.lfs),
                         climate = unique(climate),
                         continent7 = unique(continent7),
                         sphere.fragment = unique(sphere.fragment),
                         sphere.matrix = unique(sphere.matrix), 
                         biome = unique(biome),
                         taxa = unique(taxa), 
                         time.since.fragmentation = unique(time.since.fragmentation),
                         Matrix.category = unique(Matrix.category)
               ),
             by = 'dataset_label') %>% 
  # add indicator for whether study-level slope differed from zero?
  mutate(effect = ifelse((Slope_lower < 0 & Slope_upper > 0), 'random',
                         ifelse((Slope_lower > 0 & Slope_upper > 0), 'decay', 'down')))

z_Sn_studT_fragSize_group_coefs <- bind_cols(z_Sn_studT_fragSize_coef[[1]][,,'Intercept'] %>% 
                                               as_tibble() %>% 
                                               mutate(Intercept = Estimate,
                                                      Intercept_lower = Q2.5,
                                                      Intercept_upper = Q97.5,
                                                      dataset_label = rownames(z_Sn_studT_fragSize_coef[[1]][,,'Intercept'])) %>% 
                                               dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                             z_Sn_studT_fragSize_coef[[1]][,,'c.lfs'] %>% 
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
                         cxmax = max(c.lfs),
                         climate = unique(climate),
                         continent7 = unique(continent7),
                         sphere.fragment = unique(sphere.fragment),
                         sphere.matrix = unique(sphere.matrix), 
                         biome = unique(biome),
                         taxa = unique(taxa), 
                         time.since.fragmentation = unique(time.since.fragmentation),
                         Matrix.category = unique(Matrix.category)
               ),
             by = 'dataset_label') %>% 
  # add indicator for whether study-level slope differed from zero?
  mutate(effect = ifelse((Slope_lower < 0 & Slope_upper > 0), 'random',
                         ifelse((Slope_lower > 0 & Slope_upper > 0), 'decay', 'down')))

z_Scov_studT_fragSize_group_coefs <- bind_cols(z_Scov_studT_fragSize_coef[[1]][,,'Intercept'] %>% 
                                                 as_tibble() %>% 
                                                 mutate(Intercept = Estimate,
                                                        Intercept_lower = Q2.5,
                                                        Intercept_upper = Q97.5,
                                                        dataset_label = rownames(z_Scov_studT_fragSize_coef[[1]][,,'Intercept'])) %>% 
                                                 dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                               z_Scov_studT_fragSize_coef[[1]][,,'c.lfs'] %>% 
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
                         cxmax = max(c.lfs),
                         climate = unique(climate),
                         continent7 = unique(continent7),
                         sphere.fragment = unique(sphere.fragment),
                         sphere.matrix = unique(sphere.matrix), 
                         biome = unique(biome),
                         taxa = unique(taxa), 
                         time.since.fragmentation = unique(time.since.fragmentation),
                         Matrix.category = unique(Matrix.category)
               ),
             by = 'dataset_label') %>% 
  # add indicator for whether study-level slope differed from zero?
  mutate(effect = ifelse((Slope_lower < 0 & Slope_upper > 0), 'random',
                         ifelse((Slope_lower > 0 & Slope_upper > 0), 'decay', 'down')))

z_Schao_studT_fragSize_group_coefs <- bind_cols(z_Schao_studT_fragSize_coef[[1]][,,'Intercept'] %>% 
                                                  as_tibble() %>% 
                                                  mutate(Intercept = Estimate,
                                                         Intercept_lower = Q2.5,
                                                         Intercept_upper = Q97.5,
                                                         dataset_label = rownames(z_Schao_studT_fragSize_coef[[1]][,,'Intercept'])) %>% 
                                                  dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                                z_Schao_studT_fragSize_coef[[1]][,,'c.lfs'] %>% 
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
                         cxmax = max(c.lfs),
                         climate = unique(climate),
                         continent7 = unique(continent7),
                         sphere.fragment = unique(sphere.fragment),
                         sphere.matrix = unique(sphere.matrix), 
                         biome = unique(biome),
                         taxa = unique(taxa), 
                         time.since.fragmentation = unique(time.since.fragmentation),
                         Matrix.category = unique(Matrix.category)
               ),
             by = 'dataset_label') %>% 
  # add indicator for whether study-level slope differed from zero?
  mutate(effect = ifelse((Slope_lower < 0 & Slope_upper > 0), 'random',
                         ifelse((Slope_lower > 0 & Slope_upper > 0), 'decay', 'down')))

z_S_PIE_fragSize_group_coefs <- bind_cols(z_S_PIE_fS_coef[[1]][,,'Intercept'] %>% 
                                          as_tibble() %>% 
                                          mutate(Intercept = Estimate,
                                                 Intercept_lower = Q2.5,
                                                 Intercept_upper = Q97.5,
                                                 dataset_label = rownames(z_S_PIE_fS_coef[[1]][,,'Intercept'])) %>% 
                                          dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                          z_S_PIE_fS_coef[[1]][,,'c.lfs'] %>% 
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
                         cxmax = max(c.lfs),
                         climate = unique(climate),
                         continent7 = unique(continent7),
                         sphere.fragment = unique(sphere.fragment),
                         sphere.matrix = unique(sphere.matrix), 
                         biome = unique(biome),
                         taxa = unique(taxa), 
                         time.since.fragmentation = unique(time.since.fragmentation),
                         Matrix.category = unique(Matrix.category)
               ),
             by = 'dataset_label') %>% 
  # add indicator for whether study-level slope differed from zero?
  mutate(effect = ifelse((Slope_lower < 0 & Slope_upper > 0), 'random',
                         ifelse((Slope_lower > 0 & Slope_upper > 0), 'decay', 'down')))
