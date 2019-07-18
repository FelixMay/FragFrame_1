# need to execute 0_init_dirs_load_packages.R first

# code to plot visual inspection of models: trace plots, posterior predictive checks, residuals

load('~/Dropbox/1current/fragmentation_synthesis/results/jtu_z1i_fS_ref.Rdata')
load('~/Dropbox/1current/fragmentation_synthesis/results/Jne_zi_fragSize_ref.Rdata')

frag_beta <- read_csv('~/Dropbox/Frag Database (new)/files_datapaper/Analysis/2_betapart_frag_fcont_10_mabund_as_is.csv')

meta <- read.csv(paste0(path2meta, 'new_meta_2_merge.csv'), sep=';') %>% 
  as_tibble() %>% 
  dplyr::rename(dataset_label = dataset_id)

frag_beta <- left_join(frag_beta, 
                  meta,
                  by = 'dataset_label')

setwd(paste0(path2Dropbox, '/analysis_apr19/figures/model_checks'))


model = Jtu_z1i_fS
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
resid <- left_join(resid, frag_beta %>% distinct(dataset_label, sample_design),
                   by = 'dataset_label')

resid$fitted <- fitted[,'Estimate']
resid$predict <- predict[,'Estimate']
  
png('Jtu_z1i_Pearson_residuals.png', width = 240, height = 200, res = 75, units = 'mm')
par(mfrow=c(3,3), mai = c(0.5, 0.5, 0.1, 0.1))
  # can we do a better job with the ones and twos? Probably not. Error distribution? Model?
  with(resid, plot(Estimate ~ fitted,
                   ylab = 'Pearson residual'));abline(h=0, lty=2)
  with(resid, plot(Estimate ~ cl10ra,
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


