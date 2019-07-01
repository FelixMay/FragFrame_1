rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(brms)

# load data  
frag_beta <- read_csv('/gpfs1/data/idiv_chase/sablowes/fragmentation/data/sensitivity/1_betapart_frag_fcont_2_mabund_as_is.csv')
  
frag_beta <- frag_beta %>% 
    mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))

# set some weakly regularising priors
rp <- c(prior(normal(0,2), class = Intercept),
        prior(normal(0,1), class = b),
        prior(exponential(1), class = sd))
  
# fit models to baselga's components of jaccard 
Jtu_z1i_fS_1 <- brm(bf(repl ~ cl10ra + 
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

# load data  
frag_beta <- read_csv('/gpfs1/data/idiv_chase/sablowes/fragmentation/data/sensitivity/3_betapart_frag_fcont_100_mabund_as_is.csv')

frag_beta <- frag_beta %>% 
  mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))

Jtu_z1i_fS_3 <- brm(bf(repl ~ cl10ra + 
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

# load data  
frag_beta <- read_csv('/gpfs1/data/idiv_chase/sablowes/fragmentation/data/sensitivity/8_betapart_frag_fcont_10_mabund_ceiling.csv')

frag_beta <- frag_beta %>% 
  mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))

Jtu_z1i_fS_8 <- brm(bf(repl ~ cl10ra + 
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

save(Jtu_z1i_fS_1,
     Jtu_z1i_fS_3,
     Jtu_z1i_fS_8, file=Sys.getenv('OFILE'))
