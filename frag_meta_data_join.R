# check for gaps in the metadata of the fragmentation data

rm(list=ls())
library(tidyverse)

# load the latest version of the metric calculations
dat <- read_csv("~/Dropbox/Frag Database (new)/files_datapaper/Analysis/2_biodiv_frag_fcont_10_mabund_as_is.csv")

# load the old dataframe that has metadata
# meta <- read_csv('~/Dropbox/Habitat loss meta-analysis/analysis/diversity_metadata.csv')
meta <- read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
  as_tibble()

# create new column that is the first two parts of filename author_year
meta2 <- meta %>% 
  dplyr::rename(dataset_label = dataset_id)

check_merge <- left_join(dat %>% 
                    distinct(dataset_label), # we have 117 dataset_label's
                  meta2,
                  by = 'dataset_label')

check_merge %>% filter(is.na(biome)) #%>% 
  # write.csv('~/Dropbox/Frag Database (new)/files_datapaper/Analysis/missing_metadata.csv')
