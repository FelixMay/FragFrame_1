# code to plot results of fragment size (only) models
library(tidyverse)
library(brms)

load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms.Rdata')

S_std_resid <- residuals(Sstd2_lognorm_fragSize, 
                         type = 'pearson',
                         method = 'fitted') %>% 
  as_tibble() %>% 
  bind_cols(Sstd2_lognorm_fragSize$data) %>% 
  left_join(meta,
            by = 'dataset_label')

par(mfrow=c(3,3))
# can we do a better job with the ones and twos?
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
# would be good to look at the different types of data here...how what the standardisation done?