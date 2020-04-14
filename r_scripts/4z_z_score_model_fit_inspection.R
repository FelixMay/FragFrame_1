# need to execute 0_init_dirs_load_packages.R first

# code to plot visual inspection of models: trace plots, posterior predictive checks, residuals

load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_z_score_ref.Rdata')

meta <- read.csv(paste0(path2meta, 'new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag <- left_join(frag, 
                  meta,
                  by = 'dataset_label')

setwd(paste0(path2Dropbox, '/analysis_apr19/figures/model_checks'))
# create vector of response variable to loop through
response <- c('z_Sstd', 'z_S_PIE', 'z_Sn', 'z_S_chao', 'z_Scov')

for(i in 5:length(response)){
  print(paste('model', i, 'in', length(response)))
  model = paste0(response[i], '_studT_fragSize') %>% as.name()
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
  resid <- left_join(resid, frag %>% distinct(dataset_label, sample_design, z_Sstd),
                     by = 'dataset_label')
  
  resid$fitted <- fitted[,'Estimate']
  resid$predict <- predict[,'Estimate']
  
  png(paste0(response[i], '_Pearson_residuals.png'), width = 240, height = 200, res = 75, units = 'mm')
  par(mfrow=c(3,3), mai = c(0.5, 0.5, 0.1, 0.1))
  # can we do a better job with the ones and twos? Probably not. Error distribution? Model?
  with(resid,# %>% filter(Estimate<5), 
       plot(Estimate ~ fitted,
                   ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid,# %>% filter(Estimate<5), 
       plot(z_S_std ~ predict,
            ylab = 'z-score'));abline(c(0,1), lty=2)
  with(resid %>% filter(Estimate<5), 
       plot(Estimate ~ c.lfs,
                   ylab = 'Pearson residual', xlab = 'Fragment size (centred, log-scale)'));abline(h=0, lty=2)
  with(resid %>% filter(Estimate < 5),
       plot(Estimate ~ predict,
                          ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid %>% filter(Estimate<5), 
       boxplot(Estimate ~ dataset_label,
                      ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid %>% filter(Estimate<5), 
       boxplot(Estimate ~ biome,
                      ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid %>% filter(Estimate<5),
       boxplot(Estimate ~ taxa,
                      ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid %>% filter(Estimate<5), 
       boxplot(Estimate ~ Matrix.category,
                      ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid %>% filter(Estimate<5), 
       boxplot(Estimate ~ time.since.fragmentation,
                      ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid %>% filter(Estimate<5), 
       boxplot(Estimate ~ continent,
                      ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid %>% filter(Estimate<5), 
       boxplot(Estimate ~ sample_design,
                      ylab = 'Pearson residual'));abline(h=0, lty=2)
  dev.off()
}


# chain inspection
plot(SstdstudT_fragSize)
plot(NstdstudT_fragSize)
plot(S_PIEstudT_fragSize)
plot(SnstudT_fragSize)
plot(ScovstudT_fragSize)
plot(S_chaostudT_fragSize)

# posterior predictive checks
Sstd_pp <- pp_check(SstdstudT_fragSize) +
  scale_x_continuous(trans = 'log2') +
  labs(subtitle = expression(S[std]))
Nstd_pp <- pp_check(NstdstudT_fragSize) +
  scale_x_continuous(trans = 'log10') +
  labs(subtitle = expression(N[std]))
S_PIE_pp <- pp_check(S_PIEstudT_fragSize) +
  scale_x_continuous(trans = 'log2') +
  labs(subtitle = expression(S[PIE]))
Sn_pp <- pp_check(SnstudT_fragSize) +
  scale_x_continuous(trans = 'log2') +
  labs(subtitle = expression(S[n]))
Scov_pp <- pp_check(ScovstudT_fragSize) +
  scale_x_continuous(trans = 'log2') +
  labs(subtitle = expression(S[cov]))
Schao_pp <- pp_check(S_chaostudT_fragSize) +
  scale_x_continuous(trans = 'log2') +
  labs(subtitle = expression(S[chao]))

cowplot::plot_grid(Sstd_pp,
                   Nstd_pp,
                   S_PIE_pp,
                   Sn_pp,
                   Scov_pp,
                   Schao_pp,
                   nrow = 3, align = 'hv')

# ggsave('Posterior_predictive_checks.png', width = 290, height = 200, units = 'mm')
