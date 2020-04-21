# fit model to turnover components of dissimilarity (Jaccard and Ruzicka)
# execute 0_init_dirs_load_packages.R

# load the data: cluster version, then local
frag_beta <- read_csv(paste0(path2data, '2_betapart_frag_fcont_10_mabund_as_is.csv'))

# centre covariate before fitting
frag_beta <- frag_beta %>% 
  mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))


# set some weakly regularising priors
rp <- c(prior(normal(0,2), class = Intercept),
        prior(normal(0,1), class = b),
        prior(exponential(1), class = sd))

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
                        prior = rp,
                        cores = 4, chains = 4,
                        iter = 4000, thin = 2)

# save locally
save(Jtu_z1i_fragSize,
     file = '~/Dropbox/1current/fragmentation_synthesis/results/Jtu_z1i_fragSize.Rdata')

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
                        prior = rp,
                        cores = 4, chains = 4,
                        warmup = 500)

save(Rtu_z1i_fragSize,
     file = '~/Dropbox/1current/fragmentation_synthesis/results/Rtu_z1i_fragSize.Rdata')