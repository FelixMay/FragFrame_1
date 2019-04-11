# code to fit models to beta-diversity for fragmentation synthesis (fragment scale)
# code to run on rstudio server
rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(brms)

# load the data
frag_beta <- read_csv('/gpfs1/data/idiv_chase/sablowes/fragmentation/data/2_betapart_frag_fcont_10_mabund_as_is.csv')
# want to add a grouping variable for the pairwise comparisons
frag_beta <- frag_beta %>% 
  group_by(dataset_label, sample_design, method, frag_x) %>% 
  mutate(pair_group = paste0(frag_x, '_g')) %>% 
  ungroup() %>% 
  # centre covariate before fitting
  mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))


# set some weakly regularising priors
rp <- c(prior(normal(0,2), class = Intercept),
        prior(normal(0,1), class = b),
        prior(exponential(1), class = sd))

# fit models to baselga's components of jaccard 
Jtu_norm_fragSize <- brm(repl ~ cl10ra + 
                           (cl10ra | dataset_label / pair_group), 
                         # fit to data with variation in frag_size_num
                         data = frag_beta %>% filter(method=='Baselga family, Jaccard'),
                         prior = rp,
                         # family = 'lognormal',
                         cores = 4, chains = 4)


save(Jtu_norm_fragSize, file=Sys.getenv('OFILE'))