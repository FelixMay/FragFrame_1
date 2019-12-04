# execute 0_init_dirs_load_packages.R for packages and directory

# code to wrangle the coefficients for the fragemnt area regressions (with no interactions)
# for models fit to data standardised in slightly different ways

# code to get coefs for reference fits------
load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_z_score_ref.Rdata')

frag <- read_csv(paste0(path2data, '2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

meta <- read.csv(paste0(path2meta, 'new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')

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

# rename
Sstd_ref_fitted <- z_Sstd_fS_fitted
Sstd_ref_fixef <- z_Sstd_studT_fragSize_fixef
Sstd_ref_grp_coefs <- z_Sstd_studT_fragSize_group_coefs

S_PIE_ref_fitted <- z_S_PIE_fS_fitted
S_PIE_ref_fixef <- z_S_PIE_studT_fragSize_fixef
S_PIE_ref_grp_coefs <- z_S_PIE_fragSize_group_coefs


# code to get coefs for sensitivity case 1------
load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_1_biodiv_frag_fcont_2_mabund_as_is.Rdata')

frag <- read_csv(paste0(path2data, '1_biodiv_frag_fcont_2_mabund_as_is.csv'))
# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

meta <- read.csv(paste0(path2meta, 'new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')

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
# rename
Sstd_sens1_fitted <- z_Sstd_fS_fitted
Sstd_sens1_fixef <- z_Sstd_studT_fragSize_fixef
Sstd_sens1_grp_coefs <- z_Sstd_studT_fragSize_group_coefs

S_PIE_sens1_fitted <- z_S_PIE_fS_fitted
S_PIE_sens1_fixef <- z_S_PIE_studT_fragSize_fixef
S_PIE_sens1_grp_coefs <- z_S_PIE_fragSize_group_coefs


# code to get coefs for sensitivity case 3------
load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_3_biodiv_frag_fcont_100_mabund_as_is.Rdata')

frag <- read_csv(paste0(path2data, '3_biodiv_frag_fcont_100_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

meta <- read.csv(paste0(path2meta, 'new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')

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

# rename
Sstd_sens3_fitted <- z_Sstd_fS_fitted
Sstd_sens3_fixef <- z_Sstd_studT_fragSize_fixef
Sstd_sens3_grp_coefs <- z_Sstd_studT_fragSize_group_coefs

S_PIE_sens3_fitted <- z_S_PIE_fS_fitted
S_PIE_sens3_fixef <- z_S_PIE_studT_fragSize_fixef
S_PIE_sens3_grp_coefs <- z_S_PIE_fragSize_group_coefs

# code to get coefs for sensitivity case 8------
load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_4_biodiv_frag_fcont_10_mabund_ceiling.Rdata')

frag <- read_csv(paste0(path2data, '4_biodiv_frag_fcont_10_mabund_ceiling.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

meta <- read.csv(paste0(path2meta, 'new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')

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
# rename
Sstd_sens4_fitted <- z_Sstd_fS_fitted
Sstd_sens4_fixef <- z_Sstd_studT_fragSize_fixef
Sstd_sens4_grp_coefs <- z_Sstd_studT_fragSize_group_coefs

S_PIE_sens4_fitted <- z_S_PIE_fS_fitted
S_PIE_sens4_fixef <- z_S_PIE_studT_fragSize_fixef
S_PIE_sens4_grp_coefs <- z_S_PIE_fragSize_group_coefs


# code to get coefs for sensitivity case 11------
load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_5_biodiv_frag_fcont_10_mabund_multiply.Rdata')

frag <- read_csv(paste0(path2data, '5_biodiv_frag_fcont_10_mabund_multiply.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

meta <- read.csv(paste0(path2meta, 'new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')

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

# rename
Sstd_sens5_fitted <- z_Sstd_fS_fitted
Sstd_sens5_fixef <- z_Sstd_studT_fragSize_fixef
Sstd_sens5_grp_coefs <- z_Sstd_studT_fragSize_group_coefs

S_PIE_sens5_fitted <- z_S_PIE_fS_fitted
S_PIE_sens5_fixef <- z_S_PIE_studT_fragSize_fixef
S_PIE_sens5_grp_coefs <- z_S_PIE_fragSize_group_coefs