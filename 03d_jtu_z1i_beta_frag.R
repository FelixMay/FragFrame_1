rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(brms)

# execute 0_init_dirs_load_packages.R

# load the data: cluster version, then local
frag_beta <- read_csv('/gpfs1/data/idiv_chase/sablowes/fragmentation/data/2_betapart_frag_fcont_10_mabund_as_is.csv')
frag_beta <- read_csv(paste0(path2data, '2_betapart_frag_fcont_10_mabund_as_is.csv'))

# want to add a grouping variable for the pairwise comparisons
frag_beta <- frag_beta %>% 
  # centre covariate before fitting
  mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))


# set some weakly regularising priors
rp <- c(prior(normal(0,2), class = Intercept),
        prior(normal(0,1), class = b),
        prior(exponential(1), class = sd))

# fit models to baselga's components of jaccard 
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
                        cores = 4, chains = 4)

save(Jtu_z1i_fragSize, file=Sys.getenv('OFILE'))