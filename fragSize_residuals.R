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
  left_join(read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
              as_tibble() %>% 
              dplyr::rename(dataset_label = dataset_id),
            by = 'dataset_label')
  
S_std_resid4 <- residuals(Sstd2_lognorm_fragSize4, 
                         type = 'pearson',
                         method = 'fitted') %>% 
  as_tibble() %>% 
  bind_cols(Sstd2_lognorm_fragSize4$data) %>% 
  left_join(read.csv('~/Dropbox/Frag Database (new)/new_meta_2_merge.csv', sep=';') %>% 
              as_tibble() %>% 
              dplyr::rename(dataset_label = dataset_id),
            by = 'dataset_label')

# join with the sample_design column
S_std_resid <- left_join(S_std_resid, frag %>% distinct(dataset_label, sample_design),
            by = 'dataset_label')
S_std_resid4 <- left_join(S_std_resid4, frag %>% distinct(dataset_label, sample_design),
                         by = 'dataset_label')

fit1 <- fitted(Sstd2_lognorm_fragSize, re_formula = NA)
predict1 <- predict(Sstd2_lognorm_fragSize)

fit4 <- fitted(Sstd2_lognorm_fragSize4, re_formula = NA)
predict4 <- predict(Sstd2_lognorm_fragSize4)

S_std_resid$fitted <- fit1[,'Estimate']
S_std_resid$predicted <- predict1[,'Estimate']

S_std_resid4$fitted <- fit4[,'Estimate']
S_std_resid4$predicted <- predict4[,'Estimate']

par(mfrow=c(3,3))
# can we do a better job with the ones and twos? Probably not. Error distribution? Model?
with(S_std_resid, qplot(Estimate ~ climate,
                       ylab = 'Pearson residual',
                       #log = 'x'
                       ));abline(h=0, lty=2)
with(S_std_resid, plot(Estimate ~ c.lfs,
                       ylab = 'Pearson residual'));abline(h=0, lty=2)
ggplot(data = S_std_resid) +
  facet_grid(continent8 ~ time.since.fragmentation) +
  geom_boxplot(aes(y = Estimate, x = dataset_label, colour = continent8)) +
  geom_hline(yintercept = 0, lty = 2) +
  coord_flip()
       
       
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
