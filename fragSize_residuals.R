# code to plot results of fragment size (only) models
library(tidyverse)
library(brms)

load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms.Rdata')
# get the raw data (I want the sample_design column)
frag <- read_csv('~/Dropbox/Frag Database (new)/files_datapaper/Analysis/2_biodiv_frag_fcont_10_mabund_as_is.csv') %>% 
  # get the metadata 
left_join(read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id),
  by = 'dataset_label')


S_std_resid <- residuals(Sstd2_lognorm_fragSize, 
                         type = 'pearson',
                         method = 'fitted') %>% 
  as_tibble() %>% 
  bind_cols(Sstd2_lognorm_fragSize$data) %>% 
  left_join(meta,
            by = 'dataset_label')
  
# join with the sample_design column
S_std_resid <- left_join(S_std_resid, frag %>% distinct(dataset_label, sample_design),
            by = 'dataset_label')

par(mfrow=c(3,3))
# can we do a better job with the ones and twos? Probably not. Error distribution? Model?
with(S_std_resid, plot(Estimate ~ S_std_2,
                       ylab = 'Pearson residual',
                       log = 'x'));abline(h=0, lty=2)
with(S_std_resid, plot(Estimate ~ c.lfs,
                       ylab = 'Pearson residual'));abline(h=0, lty=2)
with(S_std_resid, boxplot(Estimate ~ dataset_label,
                       ylab = 'Pearson residual'));abline(h=0, lty=2)
with(S_std_resid, boxplot(Estimate ~ biome,
                          ylab = 'Pearson residual'));abline(h=0, lty=2)
with(S_std_resid, boxplot(Estimate ~ taxa,
                          ylab = 'Pearson residual'));abline(h=0, lty=2)
with(S_std_resid, boxplot(Estimate ~ Matrix.category,
                          ylab = 'Pearson residual'));abline(h=0, lty=2)
with(S_std_resid, boxplot(Estimate ~ time.since.fragmentation,
                          ylab = 'Pearson residual'));abline(h=0, lty=2)
with(S_std_resid, boxplot(Estimate ~ continent,
                          ylab = 'Pearson residual'));abline(h=0, lty=2)
with(S_std_resid, boxplot(Estimate ~ sample_design,
                          ylab = 'Pearson residual'));abline(h=0, lty=2)
