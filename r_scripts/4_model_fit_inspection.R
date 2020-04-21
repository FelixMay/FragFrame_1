# need to execute 0_init_dirs_load_packages.R first

# code to plot visual inspection of models: trace plots, posterior predictive checks, residuals
load(paste0(path2wd, '/intermediate_results/fragSize_ref.Rdata'))

meta <- read.csv(paste0(path2wd, '/data/new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag <- read_csv(paste0(path2wd, '/intermediate_results/2_biodiv_frag_fcont_10_mabund_as_is.csv'))

frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')


# change to an appropriate directory to save plots to:
plot_dir <- '~/Dropbox/1current/fragmentation_synthesis/temp/figs/visual_inspection/'
# create vector of response variable to loop through
response <- c('Sstd', 'Nstd', 'S_PIE', 'Sn', 'S_chao', 'Scov')

for(i in 5:length(response)){
  print(paste('model', i, 'in', length(response)))
  model = paste0(response[i], '_lognorm_fragSize') %>% as.name()
  resid <- residuals(eval(model), 
                          type = 'pearson',
                          method = 'predict'
                          ) %>% 
    as_tibble() %>% 
    bind_cols(eval(model)$data) %>% 
    left_join(meta,
              by = 'dataset_label')
  
  fitted <- fitted(eval(model), re_formula = NA)
  predict <- predict(eval(model))
  
  # join with the sample_design column
  resid <- left_join(resid, frag %>% distinct(dataset_label, sample_design),
                            by = 'dataset_label')
  
  resid$fitted <- fitted[,'Estimate']
  resid$predict <- predict[,'Estimate']
  
  png(paste0(plot_dir, '_', response[i], '_Pearson_residuals.png'), width = 240, height = 200, res = 75, units = 'mm')
  par(mfrow=c(3,3), mai = c(0.5, 0.5, 0.1, 0.1))
  # can we do a better job with the ones and twos? Probably not. Error distribution? Model?
  with(resid, plot(Estimate ~ fitted,
                          ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid, plot(Estimate ~ c.lfs,
                         ylab = 'Pearson residual', xlab = 'Fragment size (centred, log-scale)'));abline(h=0, lty=2)
  # with(resid, plot(Estimate ~ predict,
  #                         ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid, boxplot(Estimate ~ dataset_label,
                             ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid, boxplot(Estimate ~ biome,
                             ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid, boxplot(Estimate ~ taxa,
                             ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid, boxplot(Estimate ~ Matrix.category,
                             ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid, boxplot(Estimate ~ time.since.fragmentation,
                             ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid, boxplot(Estimate ~ continent,
                             ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid, boxplot(Estimate ~ sample_design,
                             ylab = 'Pearson residual'));abline(h=0, lty=2)
  dev.off()
}

# chain inspection
plot(Sstd_lognorm_fragSize)
plot(Nstd_lognorm_fragSize)
plot(S_PIE_lognorm_fragSize)
plot(Sn_lognorm_fragSize)
plot(Scov_lognorm_fragSize)
plot(S_chao_lognorm_fragSize)

# posterior predictive checks
Sstd_pp <- pp_check(Sstd_lognorm_fragSize) +
  scale_x_continuous(trans = 'log2') +
  labs(subtitle = expression(S[std]))
Nstd_pp <- pp_check(Nstd_lognorm_fragSize) +
  scale_x_continuous(trans = 'log10') +
  labs(subtitle = expression(N[std]))
S_PIE_pp <- pp_check(S_PIE_lognorm_fragSize) +
  scale_x_continuous(trans = 'log2') +
  labs(subtitle = expression(S[PIE]))
Sn_pp <- pp_check(Sn_lognorm_fragSize) +
  scale_x_continuous(trans = 'log2') +
  labs(subtitle = expression(S[n]))
Scov_pp <- pp_check(Scov_lognorm_fragSize) +
  scale_x_continuous(trans = 'log2') +
  labs(subtitle = expression(S[cov]))
Schao_pp <- pp_check(S_chao_lognorm_fragSize) +
  scale_x_continuous(trans = 'log2') +
  labs(subtitle = expression(S[chao]))

cowplot::plot_grid(Sstd_pp,
                   Nstd_pp,
                   S_PIE_pp,
                   Sn_pp,
                   Scov_pp,
                   Schao_pp,
                   nrow = 3, align = 'hv')

# ggsave(paste0(plot_dir, 'Posterior_predictive.png'), width = 290, height = 200, units = 'mm')
