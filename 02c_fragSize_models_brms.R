# code to fit models for fragmentation synthesis
# so far: bayesian framework for approximately ML-like results (i.e.,
# with non-informative priors - against the science of Gelman)

# 


# load the data
# frag <- read_csv('~/Dropbox/Habitat loss meta-analysis/analysis/diversity_metadata.csv')
frag <- read_csv(paste(path2temp, "diversity_metadata.csv", sep = ""))

setwd(path2temp)



# remove observations without fragment size (for now...)
frag2 <- frag %>% 
  filter(!is.na(entity.size))

# what has this filtering done to the range of entity.size?
problems <- frag2 %>% 
  group_by(filename) %>% 
  summarise(min.size = min(entity.size), 
            max.size = max(entity.size)) %>% 
  filter(min.size==max.size)
# & !is.na(taxa) & !is.na(veg.fragment) & !is.na(matrix.category) &
#            !is.na(time.since.fragmentation))

##--create some covariates for easier workflow--
# mean-centred log(fragment.size)
frag2$c.lfs <- log(frag2$entity.size) - mean(log(frag2$entity.size))
frag2$lSstd <- log(frag2$S_std)
frag2$lSn <- log(frag2$S_n)
frag2$lScov <- log(frag2$S_cov)
frag2$lS_PIE <- log(frag2$S_PIE)
frag2$lNstd <- log(frag2$N_std)

# # set the reference levels for the categorical covariates of interest
# frag2$matrix.category <- factor(frag2$matrix.category, levels = c('light filter', 'medium filter', 'harsh filter'))
# frag2$time.since.fragmentation <- factor(frag2$time.since.fragmentation,
#                                          levels = c("recent (less than 20 years)", "intermediate (20-100 years)", "long (100+ years)"))

# load model fits:
# load('~/Dropbox/1current/fragmentation_synthesis/results/fragSize_brms.Rdata')
load("fragSize_brms.Rdata")
#----- simplest model: diversity as a function of fragment size; allow fragment size to vary by study----
summary(lSstd_fragSize <- brm(lSstd ~ c.lfs + (c.lfs | filename), 
                              # fit to data with variation in entity.size
                            data = frag2 %>% filter(!filename %in% problems$filename),
                            cores = 4, chains = 4))
summary(lSstd_sknorm_fragSize <- brm(lSstd ~ c.lfs + (c.lfs | filename), 
                              family = 'skew_normal',
                              data = frag2 %>% filter(!filename %in% problems$filename),
                              cores = 4, chains = 4))
summary(lSn_fragSize <- brm(lSn ~ c.lfs + (c.lfs | filename), 
                            data = frag2 %>% filter(!filename %in% problems$filename),
                            cores = 4, chains = 4))
summary(lS_PIE_fragSize <- brm(lS_PIE ~ c.lfs + (c.lfs | filename), 
                               data = frag2 %>% filter(!filename %in% problems$filename),
                               cores = 4, chains = 4))
# note sure which abundance to model: do both
summary(N_lognorm_fragSize <- brm(N ~ c.lfs + (c.lfs | filename), 
                               data = frag2 %>% filter(!filename %in% problems$filename),
                               family = 'lognormal',
                               cores = 4, chains = 4))
summary(Nstd_lognorm_fragSize <- brm(N_std ~ c.lfs + (c.lfs | filename), 
                                  data = frag2 %>% filter(!filename %in% problems$filename &
                                                            N_std > 0),
                                  family = 'lognormal',
                                  cores = 4, chains = 4))
# refit with log-sccale response for compatibility with other models of S
summary(lNstd_fragSize <- brm(lNstd ~ c.lfs + (c.lfs | filename), 
                                     data = frag2 %>% filter(!filename %in% problems$filename &
                                                               N_std > 0),
                                     cores = 4, chains = 4))

#---quick look for consistency with data: not too bad me thinks...-----
pp_check(lSstd_fragSize)
pp_check(lSstd_sknorm_fragSize)
pp_check(lSn_fragSize)
pp_check(lS_PIE_fragSize)
pp_check(Nstd_lognorm_fragSize) +
  scale_x_continuous(trans = 'log')
pp_check(lNstd_fragSize)

save(lSstd_fragSize, lSstd_sknorm_fragSize, lSn_fragSize, lS_PIE_fragSize,
         Nstd_lognorm_fragSize, lNstd_fragSize,
      file = 'fragSize_brms.Rdata')


#------wrangle for plotting----------------
# for plotting fixed effects
lSstd_fragSize_fitted <- cbind(lSstd_fragSize$data,
                               fitted(lSstd_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag2 %>% distinct(filename, c.lfs, entity.size),
             by = c('filename', 'c.lfs'))
lSstd_fragSize_fixef <- fixef(lSstd_fragSize)

lSn_fragSize_fitted <- cbind(lSn_fragSize$data,
                               fitted(lSn_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag2 %>% distinct(filename, c.lfs, entity.size),
             by = c('filename', 'c.lfs'))

lS_PIE_fragSize_fitted <- cbind(lS_PIE_fragSize$data,
                             fitted(lS_PIE_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag2 %>% distinct(filename, c.lfs, entity.size),
             by = c('filename', 'c.lfs'))

lNstd_fragSize_fitted <- cbind(lNstd_fragSize$data,
                                fitted(lNstd_fragSize, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag2 %>% distinct(filename, c.lfs, entity.size),
             by = c('filename', 'c.lfs'))

# for plotting the random-effects
lSstd_fragSize_coef <- coef(lSstd_fragSize)
lSn_fragSize_coef <- coef(lSn_fragSize)
lS_PIE_fragSize_coef <- coef(lS_PIE_fragSize)
lNstd_fragSize_coef <- coef(lNstd_fragSize)

lSstd_fragSize_group_coefs <- bind_cols(lSstd_fragSize_coef[[1]][,,'Intercept'] %>% 
                           as_tibble() %>% 
                           mutate(Intercept = Estimate,
                                  Intercept_lower = Q2.5,
                                  Intercept_upper = Q97.5,
                                  filename = rownames(lSstd_fragSize_coef[[1]][,,'Intercept'])) %>% 
                             dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                         lSstd_fragSize_coef[[1]][,,'c.lfs'] %>% 
                           as_tibble() %>% 
                           mutate(Slope = Estimate,
                                  Slope_lower = Q2.5,
                                  Slope_upper = Q97.5) %>% 
                           dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag2 %>% 
               group_by(filename) %>% 
               summarise(xmin = min(entity.size),
                         xmax = max(entity.size),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
                         by = 'filename')

lSn_fragSize_group_coefs <- bind_cols(lSn_fragSize_coef[[1]][,,'Intercept'] %>% 
                                          as_tibble() %>% 
                                          mutate(Intercept = Estimate,
                                                 Intercept_lower = Q2.5,
                                                 Intercept_upper = Q97.5,
                                                 filename = rownames(lSn_fragSize_coef[[1]][,,'Intercept'])) %>% 
                                          dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                        lSn_fragSize_coef[[1]][,,'c.lfs'] %>% 
                                          as_tibble() %>% 
                                          mutate(Slope = Estimate,
                                                 Slope_lower = Q2.5,
                                                 Slope_upper = Q97.5) %>% 
                                          dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag2 %>% 
               group_by(filename) %>% 
               summarise(xmin = min(entity.size),
                         xmax = max(entity.size),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'filename')

lS_PIE_fragSize_group_coefs <- bind_cols(lS_PIE_fragSize_coef[[1]][,,'Intercept'] %>% 
                                          as_tibble() %>% 
                                          mutate(Intercept = Estimate,
                                                 Intercept_lower = Q2.5,
                                                 Intercept_upper = Q97.5,
                                                 filename = rownames(lS_PIE_fragSize_coef[[1]][,,'Intercept'])) %>% 
                                          dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                        lS_PIE_fragSize_coef[[1]][,,'c.lfs'] %>% 
                                          as_tibble() %>% 
                                          mutate(Slope = Estimate,
                                                 Slope_lower = Q2.5,
                                                 Slope_upper = Q97.5) %>% 
                                          dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag2 %>% 
               group_by(filename) %>% 
               summarise(xmin = min(entity.size),
                         xmax = max(entity.size),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'filename')

lNstd_fragSize_group_coefs <- bind_cols(lNstd_fragSize_coef[[1]][,,'Intercept'] %>% 
                                           as_tibble() %>% 
                                           mutate(Intercept = Estimate,
                                                  Intercept_lower = Q2.5,
                                                  Intercept_upper = Q97.5,
                                                  filename = rownames(lNstd_fragSize_coef[[1]][,,'Intercept'])) %>% 
                                           dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                       lNstd_fragSize_coef[[1]][,,'c.lfs'] %>% 
                                           as_tibble() %>% 
                                           mutate(Slope = Estimate,
                                                  Slope_lower = Q2.5,
                                                  Slope_upper = Q97.5) %>% 
                                           dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag2 %>% 
               group_by(filename) %>% 
               summarise(xmin = min(entity.size),
                         xmax = max(entity.size),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'filename')

#---- regression plots showing study-level slopes-----
# setwd('~/Dropbox/Habitat loss meta-analysis/analysis/figs/')
setwd(paste(path2temp,"figs/", sep = ""))

# S_std


ggplot() +
  # data
  geom_point(data = frag2,
             aes(x = entity.size, y = S_std, colour = filename),
             size = 1.5) +
  geom_segment(data = lSstd_fragSize_group_coefs,
               aes(group = filename,
                   colour = filename,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = lSstd_fragSize_fitted, 
               aes(x = entity.size,
                   y = exp(Estimate)),
               size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = lSstd_fragSize_fitted,
              aes(x = entity.size,
                  ymin = exp(Q2.5),
                  ymax = exp(Q97.5)),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.05, y = 132,
           label = "beta[Frag.~size] == 0.06~(0.04 - 0.09)",
           parse = T) +
  scale_x_continuous(trans = 'log', breaks = c(1e-1, 1e2, 1e5)) +
  scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[std]))) +
  theme_bw() +
  theme(legend.position = 'none')

#ggsave('S_std_regression.pdf', width = 175, height = 150, units = 'mm')
ggsave('S_std_regression.png', width = 100, height = 80, units = 'mm')


# Sn
ggplot() +
  # data
  geom_point(data = frag2,
             aes(x = entity.size, y = S_n, colour = filename),
             size = 1.5) +
  geom_segment(data = lSn_fragSize_group_coefs,
               aes(group = filename,
                   colour = filename,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = lSn_fragSize_fitted, 
            aes(x = entity.size,
                y = exp(Estimate)),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = lSn_fragSize_fitted,
              aes(x = entity.size,
                  ymin = exp(Q2.5),
                  ymax = exp(Q97.5)),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.05, y = 132,
           label = "beta[Frag.~size] == 0.05~(0.03 - 0.07)",
           parse = T) +
  scale_x_continuous(trans = 'log', breaks = c(1e-1, 1e2, 1e5)) +
  scale_y_continuous(trans = 'log', breaks = c(4,16, 32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[n]))) +
  theme_bw() +
  theme(legend.position = 'none')

# ggsave('S_n_regression.pdf', width = 175, height = 150, units = 'mm')
ggsave('S_n_regression.png', width = 100, height = 90, units = 'mm')


# S_PIE
ggplot() +
  # data
  geom_point(data = frag2,
             aes(x = entity.size, y = S_PIE, colour = filename),
             size = 1.5) +
  geom_segment(data = lS_PIE_fragSize_group_coefs,
               aes(group = filename,
                   colour = filename,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = lS_PIE_fragSize_fitted, 
            aes(x = entity.size,
                y = exp(Estimate)),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = lS_PIE_fragSize_fitted,
              aes(x = entity.size,
                  ymin = exp(Q2.5),
                  ymax = exp(Q97.5)),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  # annotate('text', x = 0.1, y = 132,
  #          label = "beta[Frag.~size] == 0.06~(0.04 - 0.08)",
  #          parse = T) +
  scale_x_continuous(trans = 'log', breaks = c(0.01,1,100,10000)) +
  scale_y_continuous(trans = 'log', breaks = c(4,16,32,64,128, 256)) +
  # scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(S[PIE]))) +
  theme_bw() +
  theme(legend.position = 'none')

# ggsave('S_PIE_regression.pdf', width = 175, height = 150, units = 'mm')
ggsave('S_PIE_regression.png', width = 100, height = 90, units = 'mm')


# N_std
ggplot() +
  # data
  geom_point(data = frag2 %>% filter(!filename %in% problems$filename &
                                                 N_std > 0),
             aes(x = entity.size, y = N_std, colour = filename),
             size = 1.5) +
  geom_segment(data = lNstd_fragSize_group_coefs,
               aes(group = filename,
                   colour = filename,
                   x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = lNstd_fragSize_fitted, 
            aes(x = entity.size,
                y = exp(Estimate)),
            size = 1.5) +
  # fixed effect uncertainty
  geom_ribbon(data = lNstd_fragSize_fitted,
              aes(x = entity.size,
                  ymin = exp(Q2.5),
                  ymax = exp(Q97.5)),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  #annotate('text', x = 0.15, y = 20000,
  #        label = "beta[Frag.~size] == 0.09~(0.03 - 0.14)",
 #parse = T, size = 4) +
  scale_x_continuous(trans = 'log10', breaks = c(0.01,1,100,10000)) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000)) +
   # scale_colour_viridis_d(guide=F) +
  labs(x = 'Fragment size (hectares)',
       y = expression(paste(N[std]))) +
  theme_bw() +
  theme(legend.position = 'none')

# ggsave('N_std_regression.pdf', width = 175, height = 150, units = 'mm')
ggsave('N_std_regression.png', width = 100, height = 90, units = 'mm')


##---coef plots---------
# combine the coefs
coefs_allRespones <- bind_rows(lS_PIE_fragSize_group_coefs %>% mutate(model = 'S_PIE'),
                               lSn_fragSize_group_coefs %>% mutate(model = 'S_n'),
                               lSstd_fragSize_group_coefs %>% mutate(model = 'S_std'),
                               Nstd_fragSize_group_coefs %>% mutate(model = 'N_std'))
# split into two for clarity
f1 <- coefs_allRespones %>% 
  distinct(filename) %>% 
  slice(1:22)
f2 <- coefs_allRespones %>% 
  distinct(filename) %>% 
  slice(23:44)
f3 <- coefs_allRespones %>% 
  distinct(filename) %>% 
  slice(45:67)
f4 <- coefs_allRespones %>% 
  distinct(filename) %>% 
  slice(68:91)

f1_panel <- ggplot() +
  # study level fragment size effects (slopes)
  geom_point(data = coefs_allRespones %>% 
               filter(filename %in% f1$filename),
             aes(x = filename,
                 y = Slope,
                 group = model,
                 colour = filename,
                 shape = model),
             position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.5) +
  geom_errorbar(data = coefs_allRespones %>% 
                  filter(filename %in% f1$filename),
                aes(x = filename,
                    ymin = Slope_lower,
                    ymax = Slope_upper,
                    group = model,
                    # linetype = model,
                    colour = filename), 
                position = position_dodge(width = 1),
                width = 0, 
                size = 0.5) +
  scale_color_viridis_d(guide = F) +
  scale_shape_manual(name = 'Response', values = c(17,19,21,12)) +
  # scale_linetype_manual(name = 'Response', values = c(1,2,3,4)) +
  # scale_x_discrete(position = 'top') +
  labs(y = expression(paste('')),
       x = '') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = c(0.8,0.1),
        panel.grid = element_blank())

f2_panel <- ggplot() +
  # study level fragment size effects (slopes)
  geom_point(data = coefs_allRespones %>% 
               filter(filename %in% f2$filename),
             aes(x = filename,
                 y = Slope,
                 group = model,
                 colour = filename,
                 shape = model),
             position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.5) +
  geom_errorbar(data = coefs_allRespones %>% 
                  filter(filename %in% f2$filename),
                aes(x = filename,
                    ymin = Slope_lower,
                    ymax = Slope_upper,
                    group = model,
                    # linetype = model,
                    colour = filename), 
                position = position_dodge(width = 1),
                width = 0, 
                size = 0.5) +
  scale_color_viridis_d(guide = F) +
  scale_shape_manual(name = 'Response', values = c(17,19,21,12)) +
  # scale_linetype_manual(name = 'Response', values = c(1,2,3,4)) +
  scale_x_discrete(position = 'top') +
  labs(y = expression(paste('')),
       x = '') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank())

f3_panel <- ggplot() +
  # study level fragment size effects (slopes)
  geom_point(data = coefs_allRespones %>% 
               filter(filename %in% f3$filename),
             aes(x = filename,
                 y = Slope,
                 group = model,
                 colour = filename,
                 shape = model),
             position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.5) +
  geom_errorbar(data = coefs_allRespones %>% 
                  filter(filename %in% f3$filename),
                aes(x = filename,
                    ymin = Slope_lower,
                    ymax = Slope_upper,
                    group = model,
                    # linetype = model,
                    colour = filename), 
                position = position_dodge(width = 1),
                width = 0, 
                size = 0.5) +
  scale_color_viridis_d(guide = F) +
  scale_shape_manual(name = 'Response', values = c(17,19,21,12)) +
  # scale_linetype_manual(name = 'Response', values = c(1,2,3,4)) +
  # scale_x_discrete(position = 'top') +
  labs(y = expression(paste('Study-level slope (fragment size)')),
       x = '') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank())

f4_panel <- ggplot() +
  # study level fragment size effects (slopes)
  geom_point(data = coefs_allRespones %>% 
               filter(filename %in% f4$filename),
             aes(x = filename,
                 y = Slope,
                 group = model,
                 colour = filename,
                 shape = model),
             position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.5) +
  geom_errorbar(data = coefs_allRespones %>% 
                  filter(filename %in% f4$filename),
                aes(x = filename,
                    ymin = Slope_lower,
                    ymax = Slope_upper,
                    group = model,
                    # linetype = model,
                    colour = filename), 
                position = position_dodge(width = 1),
                width = 0, 
                size = 0.5) +
  scale_color_viridis_d(guide = F) +
  scale_shape_manual(name = 'Response', values = c(17,19,21,12)) +
  # scale_linetype_manual(name = 'Response', values = c(1,2,3,4)) +
  scale_x_discrete(position = 'top') +
  labs(y = expression(paste('Study-level slope (fragment size)')),
       x = '') +
  coord_flip() +
  theme_bw()  +
  theme(legend.position = 'none',
        panel.grid = element_blank())

top <- cowplot::plot_grid(f1_panel, f2_panel, nrow = 1)
bottom <- cowplot::plot_grid(f3_panel, f4_panel, nrow = 1)                   
cowplot::plot_grid(top, bottom, nrow=2, align = 'hv')

# ggsave('draft_coefs.pdf', width = 300, height = 300, units = 'mm')

