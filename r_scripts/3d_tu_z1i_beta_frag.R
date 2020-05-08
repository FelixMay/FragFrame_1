# fit model to turnover components of dissimilarity (Jaccard and Ruzicka)
# execute 0_init_dirs_load_packages.R

# load the data: cluster version, then local
frag_beta <- read_csv(paste0(path2data, '2_betapart_frag_fcont_10_mabund_as_is.csv'))

# centre covariate before fitting
frag_beta <- frag_beta %>% 
  mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))


# turnover component of jaccard 
Jtu_z1i_fragSize <- brm(bf(repl ~ cl10ra + 
                          (cl10ra | dataset_label), 
                        zoi ~ cl10ra + 
                          (cl10ra | dataset_label), 
                        coi ~ cl10ra + 
                          (cl10ra | dataset_label),
                        family = zero_one_inflated_beta()),
                        # fit to data with variation in frag_size_num
                        data = frag_beta %>% filter(method=='Baselga family, Jaccard'),
                        cores = 4, chains = 4,
                        iter = 4000, thin = 2)

# save 
save(Jtu_z1i_fragSize,
     file = paste0(path2wd, 'main_results/Jtu_z1i_fragSize.Rdata'))

# turnover component of Ruzicka
Rtu_z1i_fragSize <- brm(bf(repl ~ cl10ra + 
                             (cl10ra | dataset_label), 
                           zoi ~ cl10ra + 
                             (cl10ra | dataset_label), 
                           coi ~ cl10ra + 
                             (cl10ra | dataset_label),
                           family = zero_one_inflated_beta()),
                        # fit to data with variation in frag_size_num
                        data = frag_beta %>% filter(method=='Baselga family, Ruzicka'),
                        cores = 4, chains = 4,
                        iter = 4000, thin = 2)

# save 
save(Rtu_z1i_fragSize,
     file = paste0(path2wd, 'main_results/Rtu_z1i_fragSize.Rdata'))