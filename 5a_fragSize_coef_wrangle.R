# need to execute 0_init_dirs_load_packages.R first

# code to wrangle the coefficients for the fragemnt area regressions (with no interactions)

load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms_ref_revision.Rdata')

frag <- read_csv(paste0(path2data, '2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

# load the meta data
meta <- read.csv(paste0(path2meta, 'new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')

#------wrangle for plotting
# for plotting fixed effects----------------
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
Sstd2_lognorm_fragSize_coef <- coef(Sstd2_lognorm_fragSize)
Sn_lognorm_fragSize_coef <- coef(Sn_lognorm_fragSize)
Scov_lognorm_fragSize_coef <- coef(Scov_lognorm_fragSize)
Schao_lognorm_fragSize_coef <- coef(S_chao_lognorm_fragSize)
S_PIE_fS_coef <- coef(S_PIE_lognorm_fragSize)
Nstd_fS_coef <- coef(Nstd_lognorm_fragSize)

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
