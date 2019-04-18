# code to update the models of beta-diversity to include two-way interactions
# for fragmentation synthesis (fragment scale)
# EVE
rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(brms)

# load the data
frag_beta <- read_csv('/gpfs1/data/idiv_chase/sablowes/fragmentation/data/2_betapart_frag_fcont_10_mabund_as_is.csv')

# load the meta data
meta <- read.csv('/gpfs1/data/idiv_chase/sablowes/fragmentation/data/new_meta_2_merge.csv', sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag_beta <- left_join(frag_beta, 
                       meta,
                       by = 'dataset_label')

# # to fit locally: EVE (having older version of brms package does not start sampling for the zi models of nestedness component)
# frag_beta <- read_csv('~/Dropbox/Frag Database (new)/files_datapaper/Analysis/2_betapart_frag_fcont_10_mabund_as_is.csv')
# 
# # load the meta data
# meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
#   as_tibble() %>% 
#   dplyr::rename(dataset_label = dataset_id)

# want to add a grouping variable for the pairwise comparisons
frag_beta <- frag_beta %>% 
  group_by(dataset_label, sample_design, method, frag_x) %>% 
  mutate(pair_group = paste0(frag_x, '_g')) %>% 
  ungroup() %>% 
  # centre covariate before fitting
  mutate(cl10ra = log10_ratio_area - mean(log10_ratio_area))

Rtu_z1i_fs_taxa <- brm(bf(repl ~ cl10ra*taxa + 
                            (cl10ra | dataset_label / pair_group), 
                          zoi ~ cl10ra*taxa + 
                            (cl10ra | dataset_label / pair_group), 
                          coi ~ cl10ra*taxa + 
                            (cl10ra | dataset_label / pair_group),
                          family = zero_one_inflated_beta()),
                       data = frag_beta %>% 
                         filter(method=='Baselga family, Ruzicka'),
                       cores = 4, chains = 4, iter = 2000)

save(Rtu_z1i_fs_taxa,
     file=Sys.getenv('OFILE'))
