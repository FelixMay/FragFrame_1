rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(brms)

# list of files
files = list.files(path = '/gpfs1/data/idiv_chase/sablowes/fragmentation/data/sensitivity',
                   pattern = 'betapart_frag_fcont')

frag_beta = list()

for(i in 1:length(files)){
  # get the data
  file_2_get = paste0('/gpfs1/data/idiv_chase/sablowes/fragmentation/data/sensitivity/',
                      files[i])
  
  temp <- read_csv(file_2_get)
  
  temp <- temp %>% 
    mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))
  
  frag_beta[[i]] <- temp
}

# set some weakly regularising priors
rp <- c(prior(normal(0,2), class = Intercept),
        prior(normal(0,1), class = b),
        prior(exponential(1), class = sd))

# fit models to baselga's components of Ruzicka 
Rtu_z1i_fS_1 <- brm(bf(repl ~ cl10ra + 
                         (cl10ra | dataset_label), 
                       zoi ~ cl10ra + 
                         (cl10ra | dataset_label), 
                       coi ~ cl10ra + 
                         (cl10ra | dataset_label),
                       family = zero_one_inflated_beta()),
                    # fit to data with variation in frag_size_num
                    data = frag_beta[[1]] %>% filter(method=='Baselga family, Ruzicka'),
                    prior = rp,
                    cores = 4, chains = 4)

Rtu_z1i_fS_2 <- brm(bf(repl ~ cl10ra + 
                         (cl10ra | dataset_label), 
                       zoi ~ cl10ra + 
                         (cl10ra | dataset_label), 
                       coi ~ cl10ra + 
                         (cl10ra | dataset_label),
                       family = zero_one_inflated_beta()),
                    # fit to data with variation in frag_size_num
                    data = frag_beta[[2]] %>% filter(method=='Baselga family, Ruzicka'),
                    prior = rp,
                    cores = 4, chains = 4)

Rtu_z1i_fS_3 <- brm(bf(repl ~ cl10ra + 
                         (cl10ra | dataset_label), 
                       zoi ~ cl10ra + 
                         (cl10ra | dataset_label), 
                       coi ~ cl10ra + 
                         (cl10ra | dataset_label),
                       family = zero_one_inflated_beta()),
                    # fit to data with variation in frag_size_num
                    data = frag_beta[[3]] %>% filter(method=='Baselga family, Ruzicka'),
                    prior = rp,
                    cores = 4, chains = 4)

Rtu_z1i_fS_4 <- brm(bf(repl ~ cl10ra + 
                         (cl10ra | dataset_label), 
                       zoi ~ cl10ra + 
                         (cl10ra | dataset_label), 
                       coi ~ cl10ra + 
                         (cl10ra | dataset_label),
                       family = zero_one_inflated_beta()),
                    # fit to data with variation in frag_size_num
                    data = frag_beta[[4]] %>% filter(method=='Baselga family, Ruzicka'),
                    prior = rp,
                    cores = 4, chains = 4)

names = c(strsplit(files[1], split = '.csv')[[1]],
          strsplit(files[2], split = '.csv')[[1]],
          strsplit(files[3], split = '.csv')[[1]],
          strsplit(files[4], split = '.csv')[[1]])

save(names, 
     Rtu_z1i_fS_1,
     Rtu_z1i_fS_2,
     Rtu_z1i_fS_3,
     Rtu_z1i_fS_4, file=Sys.getenv('OFILE'))
