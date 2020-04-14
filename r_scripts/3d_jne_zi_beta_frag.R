# code to fit models to beta-diversity for fragmentation synthesis (fragment scale)


# code to run on rstudio server or EVE 
rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(brms)

# alt, if running locally execute 0_init_dirs_load_packages.R
# load the data: cluster version
frag_beta <- read_csv('/gpfs1/data/idiv_chase/sablowes/fragmentation/data/2_betapart_frag_fcont_10_mabund_as_is.csv')
# to fit locally: EVE (having older version of brms package does not start sampling for the zi models of nestedness component)
frag_beta <- read_csv(paste0(path2data, '2_betapart_frag_fcont_10_mabund_as_is.csv'))

# centre covariate before fitting
frag_beta <- frag_beta %>% 
  mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))


# set some weakly regularising priors
get_prior(bf(rich ~ cl10ra + 
               (cl10ra | dataset_label / pair_group), 
             zi ~ cl10ra + 
               (cl10ra | dataset_label / pair_group),
             family = zero_inflated_beta()),
          data = frag_beta %>% filter(method=='Baselga family, Jaccard'))

rp <- c(prior(normal(0,2), class = Intercept),
        prior(logistic(0,1), class = Intercept, dpar = zi),
        prior(normal(0,1), class = b),
        prior(exponential(1), class = sd))

# fit models to baselga's components of jaccard 
Jne_zi_fragSize <- brm(bf(rich ~ cl10ra + 
                             (cl10ra | dataset_label), 
                           zi ~ cl10ra + 
                             (cl10ra | dataset_label),
                           family = zero_inflated_beta()),
                        # fit to data with variation in frag_size_num
                        data = frag_beta %>% filter(method=='Baselga family, Jaccard'),
                        prior = rp,
                        cores = 4, chains = 4, iter = 2000)

save(Jne_zi_fragSize, file='~/Dropbox/1current/fragmentation_synthesis/results/Jne_zi_fragSize_ref.Rdata')
# plot(Jne_zi_fragSize)
