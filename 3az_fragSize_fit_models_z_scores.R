# new code for revision: fit models to z-scores (instead of the raw diversity measures) 
# a z-score for each diversity response (not the numbers of individuals, N) was calculated as:
# z = (observed - expected) / sd(expected)

# for this analysis I am sticking with the default priors for simplicity 
# they are weakly regularising

# load the data
frag <- read_csv(paste0(path2data, '2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

# visual inspection of the z-scores and posterior predictive checks of models fit with gaussian
# error show a big spike in z-scores around zero. Best chance of modelling these (I think) is with
# the asymmetric laplace distribution or student T


# new addition for revision: model z-scores: (observed(S_std) - expected(S_Std)/sd(expected(S_std)))
z_Sstd_asymL_fragSize <- brm(z_S_std ~ c.lfs + (c.lfs | dataset_label), 
                             # some z-scores are infinite due to sd(expected) = 0
                             data = frag %>% filter(S_std_mean>0 & 
                                                      !is.na(z_S_std) & 
                                                      !is.infinite(z_S_std)),
                             # prior = rp,
                             # our z-scores are distributed around zero so leave family as default gaussian
                             # alternate use asymmetric laplace to do better job of modelling spike in z-scores
                             # centred at zero
                             family = asym_laplace(),
                             cores = 4, 
                             chains = 4
)

z_Sstd_studT_fragSize <- brm(z_S_std ~ c.lfs + (c.lfs | dataset_label), 
                             # some z-scores are infinite due to sd(expected) = 0
                             data = frag %>% filter(S_std_mean>0 & 
                                                      !is.na(z_S_std) & 
                                                      !is.infinite(z_S_std)),
                             # prior = rp,
                             # our z-scores are distributed around zero so leave family as default gaussian
                             # alternate use asymmetric laplace to do better job of modelling spike in z-scores
                             # centred at zero
                             family = student(),
                             cores = 4, 
                             chains = 4
)


z_Sstd_studT_fragSize <- add_criterion(z_Sstd_studT_fragSize, criterion = 'loo')
z_Sstd_asymL_fragSize <- add_criterion(z_Sstd_asymL_fragSize, criterion = 'loo')

loo::loo_compare(z_Sstd_studT_fragSize,
            z_Sstd_asymL_fragSize,
            criterion = 'loo')

model_weights(z_Sstd_studT_fragSize,
            z_Sstd_asymL_fragSize,
            weights = 'loo')

z_Sn_asymL_fragSize <- brm(z_S_n ~ c.lfs + (c.lfs | dataset_label), 
                           # some z-scores are infinite due to sd(expected) = 0
                           data = frag %>% filter(S_n_mean>0 & 
                                                    !is.na(z_S_n) & 
                                                    !is.infinite(z_S_n)),
                           # prior = rp,
                           # our z-scores are distributed around zero so leave family as default gaussian
                           # alternate use asymmetric laplace to do better job of modelling spike in z-scores
                           # centred at zero
                           family = asym_laplace(),
                           cores = 4, 
                           chains = 4
)

z_Sn_studT_fragSize <- brm(z_S_n ~ c.lfs + (c.lfs | dataset_label), 
                           # some z-scores are infinite due to sd(expected) = 0
                           data = frag %>% filter(S_n_mean>0 & 
                                                    !is.na(z_S_n) & 
                                                    !is.infinite(z_S_n)),
                           # prior = rp,
                           # spike in values near zero is not well described with gaussian error
                           # alternate use student T
                           # centred at zero
                           family = student(),
                           cores = 4, 
                           chains = 4
)

z_Sn_skewN_fragSize <- brm(z_S_n ~ c.lfs + (c.lfs | dataset_label), 
                           # some z-scores are infinite due to sd(expected) = 0
                           data = frag %>% filter(S_n_mean>0 & 
                                                    !is.na(z_S_n) & 
                                                    !is.infinite(z_S_n)),
                           # prior = rp,
                           # spike in values near zero is not well described with gaussian error
                           # alternate use student T
                           family = skew_normal(),
                           cores = 4, 
                           chains = 4
)

pp_check(z_Sn_asymL_fragSize) +
  scale_x_continuous(limits = c(-10, 10))

pp_check(z_Sn_studT_fragSize) +
  scale_x_continuous(limits = c(-10, 10))

pp_check(z_Sn_skewN_fragSize) +
  scale_x_continuous(limits = c(-10, 10))

z_Sn_asymL_fragSize <- add_criterion(z_Sn_asymL_fragSize, criterion = 'loo')
z_Sn_studT_fragSize <- add_criterion(z_Sn_studT_fragSize, criterion = 'loo')
z_Sn_skewN_fragSize <- add_criterion(z_Sn_skewN_fragSize, criterion = 'loo')

loo::loo_compare(z_Sn_asymL_fragSize,
                 z_Sn_studT_fragSize,
                 z_Sn_skewN_fragSize,
                 criterion = 'loo')

model_weights(z_Sn_studT_fragSize,
              z_Sn_asymL_fragSize,
              z_Sn_skewN_fragSize,
              weights = 'loo')

z_Scov_asymL_fragSize <- brm(z_S_cov ~ c.lfs + (c.lfs | dataset_label), 
                             # some z-scores are infinite due to sd(expected) = 0
                             data = frag %>% filter(S_cov_mean>0 & 
                                                      !is.na(z_S_cov) & 
                                                      !is.infinite(z_S_cov)),
                             # prior = rp,
                             # our z-scores are distributed around zero so leave family as default gaussian
                             # alternate use asymmetric laplace to do better job of modelling spike in z-scores
                             # centred at zero
                             family = asym_laplace(),
                             cores = 4, 
                             chains = 4
)

z_Scov_studT_fragSize <- brm(z_S_cov ~ c.lfs + (c.lfs | dataset_label), 
                             # some z-scores are infinite due to sd(expected) = 0
                             data = frag %>% filter(S_cov_mean>0 & 
                                                      !is.na(z_S_cov) & 
                                                      !is.infinite(z_S_cov)),
                             # prior = rp,
                             # spike in values near zero is not well described with gaussian error
                             # alternate use student T
                             family = student(),
                             cores = 4, 
                             chains = 4
)

pp_check(z_Scov_asymL_fragSize) +
  scale_x_continuous(limits = c(-10, 10))

pp_check(z_Scov_studT_fragSize) +
  scale_x_continuous(limits = c(-10, 10))

z_Scov_asymL_fragSize <- add_criterion(z_Scov_asymL_fragSize, criterion = 'loo')
z_Scov_studT_fragSize <- add_criterion(z_Scov_studT_fragSize, criterion = 'loo')


loo::loo_compare(z_Scov_asymL_fragSize,
                 z_Scov_studT_fragSize,
                 criterion = 'loo')

model_weights(z_Scov_asymL_fragSize,
              z_Scov_studT_fragSize,
              weights = 'loo')

# S_PIE
z_S_PIE_asymL_fragSize <- brm(z_S_PIE ~ c.lfs + (c.lfs | dataset_label), 
                             # some z-scores are infinite due to sd(expected) = 0
                             data = frag %>% filter(S_PIE_mean>0 & 
                                                      !is.na(z_S_PIE) & 
                                                      !is.infinite(z_S_PIE)),
                             # prior = rp,
                             # spike in values near zero is not well described with gaussian error
                             # alternate use asymmetric laplace
                             family = asym_laplace(),
                             cores = 4, 
                             chains = 4
)

z_S_PIE_studT_fragSize <- brm(z_S_PIE ~ c.lfs + (c.lfs | dataset_label), 
                              # some z-scores are infinite due to sd(expected) = 0
                              data = frag %>% filter(S_PIE_mean>0 & 
                                                       !is.na(z_S_PIE) & 
                                                       !is.infinite(z_S_PIE)),
                              # prior = rp,
                              # spike in values near zero is not well described with gaussian error
                              # alternate use student T
                              family = student(),
                              cores = 4, 
                              chains = 4
)


pp_check(z_S_PIE_asymL_fragSize) +
  scale_x_continuous(limits = c(-10, 10))

pp_check(z_S_PIE_studT_fragSize) +
  scale_x_continuous(limits = c(-10, 10))

z_S_PIE_asymL_fragSize <- add_criterion(z_S_PIE_asymL_fragSize, criterion = 'loo')
z_S_PIE_studT_fragSize <- add_criterion(z_S_PIE_studT_fragSize, criterion = 'loo')


loo::loo_compare(z_S_PIE_asymL_fragSize,
                 z_S_PIE_studT_fragSize,
                 criterion = 'loo')

model_weights(z_S_PIE_asymL_fragSize,
              z_S_PIE_studT_fragSize,
              weights = 'loo')

# S_chao
z_S_chao_asymL_fragSize <- brm(z_S_chao ~ c.lfs + (c.lfs | dataset_label), 
                              # some z-scores are infinite due to sd(expected) = 0
                              data = frag %>% filter(S_chao_mean>0 & 
                                                       !is.na(z_S_chao) & 
                                                       !is.infinite(z_S_chao)),
                              # prior = rp,
                              # spike in values near zero is not well described with gaussian error
                              # alternate use asymmetric laplace
                              family = asym_laplace(),
                              cores = 4, 
                              chains = 4
)

z_S_chao_studT_fragSize <- brm(z_S_chao ~ c.lfs + (c.lfs | dataset_label), 
                              # some z-scores are infinite due to sd(expected) = 0
                              data = frag %>% filter(S_chao_mean>0 & 
                                                       !is.na(z_S_chao) & 
                                                       !is.infinite(z_S_chao)),
                              # prior = rp,
                              # spike in values near zero is not well described with gaussian error
                              # alternate use student T
                              family = student(),
                              cores = 4, 
                              chains = 4
)


pp_check(z_S_chao_asymL_fragSize) +
  scale_x_continuous(limits = c(-10, 10))

pp_check(z_S_chao_studT_fragSize) +
  scale_x_continuous(limits = c(-10, 10))

z_S_chao_asymL_fragSize <- add_criterion(z_S_chao_asymL_fragSize, criterion = 'loo')
z_S_chao_studT_fragSize <- add_criterion(z_S_chao_studT_fragSize, criterion = 'loo')


loo::loo_compare(z_S_chao_asymL_fragSize,
                 z_S_chao_studT_fragSize,
                 criterion = 'loo')

model_weights(z_S_chao_asymL_fragSize,
              z_S_chao_studT_fragSize,
              weights = 'loo')


save(z_Sstd_studT_fragSize,
     z_Sn_studT_fragSize,
     z_Scov_studT_fragSize,
     z_S_PIE_studT_fragSize,
     z_S_chao_studT_fragSize,
     file = '~/Dropbox/1current/fragmentation_synthesis/results/fragSize_z_score_ref.Rdata')
