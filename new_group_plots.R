# code to look at the cbntt coefficients

Sstd2_fS_fitted4 <- cbind(Sstd2_lognorm_fragSize4$data,
                         fitted(Sstd2_lognorm_fragSize4, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(cbmtt, c.lfs, frag_size_num),
             by = c('cbmtt', 'c.lfs'))

Sstd2_lognorm_fragSize_fixef4 <- fixef(Sstd2_lognorm_fragSize4)

Sn_fS_fitted4 <- cbind(Sn_lognorm_fragSize4$data,
                      fitted(Sn_lognorm_fragSize4, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(cbmtt, c.lfs, frag_size_num),
             by = c('cbmtt', 'c.lfs'))

Sn_lognorm_fragSize_fixef4 <- fixef(Sn_lognorm_fragSize4)

Scov_fS_fitted4 <- cbind(Scov_lognorm_fragSize4$data,
                        fitted(Scov_lognorm_fragSize4, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(cbmtt, c.lfs, frag_size_num),
             by = c('cbmtt', 'c.lfs'))

Scov_lognorm_fragSize_fixef4 <- fixef(Scov_lognorm_fragSize4)

Schao_fS_fitted4 <- cbind(S_chao_lognorm_fragSize4$data,
                         fitted(S_chao_lognorm_fragSize4, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(cbmtt, c.lfs, frag_size_num),
             by = c('cbmtt', 'c.lfs'))

Schao_lognorm_fragSize_fixef4 <- fixef(S_chao_lognorm_fragSize4)

S_PIE_fS_fitted4 <- cbind(S_PIE_lognorm_fragSize4$data,
                         fitted(S_PIE_lognorm_fragSize4, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(cbmtt, c.lfs, frag_size_num),
             by = c('cbmtt', 'c.lfs'))

S_PIE_lognorm_fragSize_fixef4 <- fixef(S_PIE_lognorm_fragSize4)

Nstd_fS_fitted4 <- cbind(Nstd_lognorm_fragSize4$data,
                        fitted(Nstd_lognorm_fragSize4, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(frag %>% distinct(cbmtt, c.lfs, frag_size_num),
             by = c('cbmtt', 'c.lfs'))

Nstd_lognorm_fragSize_fixef4 <- fixef(Nstd_lognorm_fragSize4)


# for plotting the random-effects----------------
Sstd2_lognorm_fragSize_coef4 <- coef(Sstd2_lognorm_fragSize4)
Sn_lognorm_fragSize_coef4 <- coef(Sn_lognorm_fragSize4)
Scov_lognorm_fragSize_coef4 <- coef(Scov_lognorm_fragSize4)
Schao_lognorm_fragSize_coef4 <- coef(S_chao_lognorm_fragSize4)
S_PIE_fS_coef4 <- coef(S_PIE_lognorm_fragSize4)
Nstd_fS_coef4 <- coef(Nstd_lognorm_fragSize4)

Sstd2_lognorm_fragSize_group_coefs4 <- bind_cols(Sstd2_lognorm_fragSize_coef4[[1]][,,'Intercept'] %>% 
                                                  as_tibble() %>% 
                                                  mutate(Intercept = Estimate,
                                                         Intercept_lower = Q2.5,
                                                         Intercept_upper = Q97.5,
                                                         cbmtt = rownames(Sstd2_lognorm_fragSize_coef4[[1]][,,'Intercept'])) %>% 
                                                  dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                                Sstd2_lognorm_fragSize_coef4[[1]][,,'c.lfs'] %>% 
                                                  as_tibble() %>% 
                                                  mutate(Slope = Estimate,
                                                         Slope_lower = Q2.5,
                                                         Slope_upper = Q97.5) %>% 
                                                  dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(cbmtt) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'cbmtt')

Sn_lognorm_fragSize_group_coefs4 <- bind_cols(Sn_lognorm_fragSize_coef4[[1]][,,'Intercept'] %>% 
                                               as_tibble() %>% 
                                               mutate(Intercept = Estimate,
                                                      Intercept_lower = Q2.5,
                                                      Intercept_upper = Q97.5,
                                                      cbmtt = rownames(Sn_lognorm_fragSize_coef4[[1]][,,'Intercept'])) %>% 
                                               dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                             Sn_lognorm_fragSize_coef4[[1]][,,'c.lfs'] %>% 
                                               as_tibble() %>% 
                                               mutate(Slope = Estimate,
                                                      Slope_lower = Q2.5,
                                                      Slope_upper = Q97.5) %>% 
                                               dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(cbmtt) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'cbmtt')

Scov_lognorm_fragSize_group_coefs4 <- bind_cols(Scov_lognorm_fragSize_coef4[[1]][,,'Intercept'] %>% 
                                                 as_tibble() %>% 
                                                 mutate(Intercept = Estimate,
                                                        Intercept_lower = Q2.5,
                                                        Intercept_upper = Q97.5,
                                                        cbmtt = rownames(Scov_lognorm_fragSize_coef4[[1]][,,'Intercept'])) %>% 
                                                 dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                               Scov_lognorm_fragSize_coef4[[1]][,,'c.lfs'] %>% 
                                                 as_tibble() %>% 
                                                 mutate(Slope = Estimate,
                                                        Slope_lower = Q2.5,
                                                        Slope_upper = Q97.5) %>% 
                                                 dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(cbmtt) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'cbmtt')

Schao_lognorm_fragSize_group_coefs4 <- bind_cols(Schao_lognorm_fragSize_coef4[[1]][,,'Intercept'] %>% 
                                                  as_tibble() %>% 
                                                  mutate(Intercept = Estimate,
                                                         Intercept_lower = Q2.5,
                                                         Intercept_upper = Q97.5,
                                                         cbmtt = rownames(Schao_lognorm_fragSize_coef4[[1]][,,'Intercept'])) %>% 
                                                  dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                                Schao_lognorm_fragSize_coef4[[1]][,,'c.lfs'] %>% 
                                                  as_tibble() %>% 
                                                  mutate(Slope = Estimate,
                                                         Slope_lower = Q2.5,
                                                         Slope_upper = Q97.5) %>% 
                                                  dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(cbmtt) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'cbmtt')

S_PIE_fragSize_group_coefs4 <- bind_cols(S_PIE_fS_coef4[[1]][,,'Intercept'] %>% 
                                          as_tibble() %>% 
                                          mutate(Intercept = Estimate,
                                                 Intercept_lower = Q2.5,
                                                 Intercept_upper = Q97.5,
                                                 cbmtt = rownames(S_PIE_fS_coef4[[1]][,,'Intercept'])) %>% 
                                          dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                        S_PIE_fS_coef4[[1]][,,'c.lfs'] %>% 
                                          as_tibble() %>% 
                                          mutate(Slope = Estimate,
                                                 Slope_lower = Q2.5,
                                                 Slope_upper = Q97.5) %>% 
                                          dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(cbmtt) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'cbmtt')

Nstd_fragSize_group_coefs4 <- bind_cols(Nstd_fS_coef4[[1]][,,'Intercept'] %>% 
                                         as_tibble() %>% 
                                         mutate(Intercept = Estimate,
                                                Intercept_lower = Q2.5,
                                                Intercept_upper = Q97.5,
                                                cbmtt = rownames(Nstd_fS_coef4[[1]][,,'Intercept'])) %>% 
                                         dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                       Nstd_fS_coef4[[1]][,,'c.lfs'] %>% 
                                         as_tibble() %>% 
                                         mutate(Slope = Estimate,
                                                Slope_lower = Q2.5,
                                                Slope_upper = Q97.5) %>% 
                                         dplyr::select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(frag %>% 
               group_by(cbmtt) %>% 
               summarise(xmin = min(frag_size_num),
                         xmax = max(frag_size_num),
                         cxmin = min(c.lfs),
                         cxmax = max(c.lfs)),
             by = 'cbmtt')


Sstd2_lognorm_fragSize_group_coefs4 <- Sstd2_lognorm_fragSize_group_coefs4 %>% 
  separate(col = cbmtt, into = c('Continent', 'Biome', 'Matrix', 'Time', 'Taxa'), remove = F, sep = '_') 
Sstd2_lognorm_fragSize_group_coefs4$Matrix <- factor(Sstd2_lognorm_fragSize_group_coefs4$Matrix,
                                                     levels = c('light filter', 'intermediate', 'harsh filter'))
Sstd2_lognorm_fragSize_group_coefs4$Time <- factor(Sstd2_lognorm_fragSize_group_coefs4$Time,
                                                     levels = c('Recent (less than 20 years)', 'Intermediate (20-100 years)', 'long (100+ years)'),
                                                  labels = c('< 20 years', '20-100 years', '> 100 years'))

Scov_lognorm_fragSize_group_coefs4 <- Scov_lognorm_fragSize_group_coefs4 %>% 
  separate(col = cbmtt, into = c('Continent', 'Biome', 'Matrix', 'Time', 'Taxa'), remove = F, sep = '_') 
Scov_lognorm_fragSize_group_coefs4$Matrix <- factor(Scov_lognorm_fragSize_group_coefs4$Matrix,
                                                     levels = c('light filter', 'intermediate', 'harsh filter'))
Scov_lognorm_fragSize_group_coefs4$Time <- factor(Scov_lognorm_fragSize_group_coefs4$Time,
                                                   levels = c('Recent (less than 20 years)', 'Intermediate (20-100 years)', 'long (100+ years)'),
                                                   labels = c('< 20 years', '20-100 years', '> 100 years'))

Sn_lognorm_fragSize_group_coefs4 <- Sn_lognorm_fragSize_group_coefs4 %>% 
  separate(col = cbmtt, into = c('Continent', 'Biome', 'Matrix', 'Time', 'Taxa'), remove = F, sep = '_') 
Sn_lognorm_fragSize_group_coefs4$Matrix <- factor(Sn_lognorm_fragSize_group_coefs4$Matrix,
                                                     levels = c('light filter', 'intermediate', 'harsh filter'))
Sn_lognorm_fragSize_group_coefs4$Time <- factor(Sn_lognorm_fragSize_group_coefs4$Time,
                                                   levels = c('Recent (less than 20 years)', 'Intermediate (20-100 years)', 'long (100+ years)'),
                                                   labels = c('< 20 years', '20-100 years', '> 100 years'))

Schao_lognorm_fragSize_group_coefs4 <- Schao_lognorm_fragSize_group_coefs4 %>% 
  separate(col = cbmtt, into = c('Continent', 'Biome', 'Matrix', 'Time', 'Taxa'), remove = F, sep = '_') 
Schao_lognorm_fragSize_group_coefs4$Matrix <- factor(Schao_lognorm_fragSize_group_coefs4$Matrix,
                                                     levels = c('light filter', 'intermediate', 'harsh filter'))

Schao_lognorm_fragSize_group_coefs4$Time <- factor(Schao_lognorm_fragSize_group_coefs4$Time,
                                                   levels = c('Recent (less than 20 years)', 'Intermediate (20-100 years)', 'long (100+ years)'),
                                                   labels = c('< 20 years', '20-100 years', '> 100 years'))

S_PIE_fragSize_group_coefs4 <- S_PIE_fragSize_group_coefs4 %>% 
  separate(col = cbmtt, into = c('Continent', 'Biome', 'Matrix', 'Time', 'Taxa'), remove = F, sep = '_') 
S_PIE_fragSize_group_coefs4$Matrix <- factor(S_PIE_fragSize_group_coefs4$Matrix,
                                                     levels = c('light filter', 'intermediate', 'harsh filter'))
S_PIE_fragSize_group_coefs4$Time <- factor(S_PIE_fragSize_group_coefs4$Time,
                                                   levels = c('Recent (less than 20 years)', 'Intermediate (20-100 years)', 'long (100+ years)'),
                                                   labels = c('< 20 years', '20-100 years', '> 100 years'))

Nstd_fragSize_group_coefs4 <- Nstd_fragSize_group_coefs4 %>% 
  separate(col = cbmtt, into = c('Continent', 'Biome', 'Matrix', 'Time', 'Taxa'), remove = F, sep = '_') 
Nstd_fragSize_group_coefs4$Matrix <- factor(Nstd_fragSize_group_coefs4$Matrix,
                                                     levels = c('light filter', 'intermediate', 'harsh filter'))
Nstd_fragSize_group_coefs4$Time <- factor(Nstd_fragSize_group_coefs4$Time,
                                                   levels = c('Recent (less than 20 years)', 'Intermediate (20-100 years)', 'long (100+ years)'),
                                                   labels = c('< 20 years', '20-100 years', '> 100 years'))


ggplot() +
  facet_grid(Continent ~ Taxa) +
  geom_point(data = Sstd2_lognorm_fragSize_group_coefs4,
             aes(x = Biome, y = Slope, colour = Matrix, shape = Time),
             position = position_dodge(width = 1),
             size = 2) +
  geom_linerange(data = Sstd2_lognorm_fragSize_group_coefs4, 
                 aes(x = Biome,
                     ymin = Slope_lower, ymax = Slope_upper, 
                     colour = Matrix, linetype = Time),
                 position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0, size = 0.5, lty = 2) +
  geom_hline(data = as.data.frame(Sstd2_lognorm_fragSize_fixef4),
             aes(yintercept = Estimate[2])) +
  geom_rect(data = as.data.frame(Sstd2_lognorm_fragSize_fixef4),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(subtitle = expression(paste(S[std]))) +
  coord_flip() +
  theme_bw()
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/Sstd_new_group.pdf', width = 290, height = 200, units = 'mm')

ggplot() +
  facet_grid(Continent ~ Taxa) +
  geom_point(data = Scov_lognorm_fragSize_group_coefs4,
             aes(x = Biome, y = Slope, colour = Matrix, shape = Time),
             position = position_dodge(width = 1),
             size = 2) +
  geom_linerange(data = Scov_lognorm_fragSize_group_coefs4, 
                 aes(x = Biome,
                     ymin = Slope_lower, ymax = Slope_upper, 
                     colour = Matrix, linetype = Time),
                 position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0, size = 0.5, lty = 2) +
  geom_hline(data = as.data.frame(Scov_lognorm_fragSize_fixef4),
             aes(yintercept = Estimate[2])) +
  geom_rect(data = as.data.frame(Scov_lognorm_fragSize_fixef4),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(subtitle = expression(paste(S[cov]))) +
  coord_flip() +
  theme_bw()
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/Scov_new_group.pdf', width = 290, height = 200, units = 'mm')

ggplot() +
  facet_grid(Continent ~ Taxa) +
  geom_point(data = Sn_lognorm_fragSize_group_coefs4,
             aes(x = Biome, y = Slope, colour = Matrix, shape = Time),
             position = position_dodge(width = 1),
             size = 2) +
  geom_linerange(data = Sn_lognorm_fragSize_group_coefs4, 
                 aes(x = Biome,
                     ymin = Slope_lower, ymax = Slope_upper, 
                     colour = Matrix, linetype = Time),
                 position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0, size = 0.5, lty = 2) +
  geom_hline(data = as.data.frame(Sn_lognorm_fragSize_fixef4),
             aes(yintercept = Estimate[2])) +
  geom_rect(data = as.data.frame(Sn_lognorm_fragSize_fixef4),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(subtitle = expression(paste(S[n]))) +
  coord_flip() +
  theme_bw()
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/Sn_new_group.pdf', width = 290, height = 200, units = 'mm')

ggplot() +
  facet_grid(Continent ~ Taxa) +
  geom_point(data = Schao_lognorm_fragSize_group_coefs4,
             aes(x = Biome, y = Slope, colour = Matrix, shape = Time),
             position = position_dodge(width = 1),
             size = 2) +
  geom_linerange(data = Schao_lognorm_fragSize_group_coefs4, 
                 aes(x = Biome,
                     ymin = Slope_lower, ymax = Slope_upper, 
                     colour = Matrix, linetype = Time),
                 position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0, size = 0.5, lty = 2) +
  geom_hline(data = as.data.frame(Schao_lognorm_fragSize_fixef4),
             aes(yintercept = Estimate[2])) +
  geom_rect(data = as.data.frame(Schao_lognorm_fragSize_fixef4),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(subtitle = expression(paste(S[chao]))) +
  coord_flip() +
  theme_bw()
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/Schao_new_group.pdf', width = 290, height = 200, units = 'mm')

ggplot() +
  facet_grid(Continent ~ Taxa) +
  geom_point(data = S_PIE_fragSize_group_coefs4,
             aes(x = Biome, y = Slope, colour = Matrix, shape = Time),
             position = position_dodge(width = 1),
             size = 2) +
  geom_linerange(data = S_PIE_fragSize_group_coefs4, 
                 aes(x = Biome,
                     ymin = Slope_lower, ymax = Slope_upper, 
                     colour = Matrix, linetype = Time),
                 position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0, size = 0.5, lty = 2) +
  geom_hline(data = as.data.frame(S_PIE_lognorm_fragSize_fixef4),
             aes(yintercept = Estimate[2])) +
  geom_rect(data = as.data.frame(S_PIE_lognorm_fragSize_fixef4),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(subtitle = expression(paste(S[PIE]))) +
  coord_flip() +
  theme_bw()
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/Spie_new_group.pdf', width = 290, height = 200, units = 'mm')

ggplot() +
  facet_grid(Continent ~ Taxa) +
  geom_point(data = Nstd_fragSize_group_coefs4,
             aes(x = Biome, y = Slope, colour = Matrix, shape = Time),
             position = position_dodge(width = 1),
             size = 2) +
  geom_linerange(data = Nstd_fragSize_group_coefs4, 
                 aes(x = Biome,
                     ymin = Slope_lower, ymax = Slope_upper, 
                     colour = Matrix, linetype = Time),
                 position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0, size = 0.5, lty = 2) +
  geom_hline(data = as.data.frame(Nstd_lognorm_fragSize_fixef4),
             aes(yintercept = Estimate[2])) +
  geom_rect(data = as.data.frame(Nstd_lognorm_fragSize_fixef4),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(subtitle = expression(paste(N[std]))) +
  coord_flip() +
  theme_bw()
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/Nstd_new_group.pdf', width = 290, height = 200, units = 'mm')


ggplot() +
  facet_grid(Continent ~ Taxa) +
  geom_point(data = Sstd2_lognorm_fragSize_group_coefs4 %>% filter(Biome=='forest'),
             aes(x = Time, y = Slope, colour = Matrix, group = interaction(Matrix,Continent)),
             position = position_dodge(width = 0.5),
             size = 2) +
  geom_linerange(data = Sstd2_lognorm_fragSize_group_coefs4 %>% filter(Biome=='forest'), 
                 aes(x = Time,
                     ymin = Slope_lower, ymax = Slope_upper, 
                     colour = Matrix, group = interaction(Matrix,Continent)),
                 position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, size = 0.5, lty = 2) +
  geom_hline(data = as.data.frame(Sstd2_lognorm_fragSize_fixef4),
             aes(yintercept = Estimate[2])) +
  geom_rect(data = as.data.frame(Sstd2_lognorm_fragSize_fixef4),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(subtitle = expression(paste(S[std]))) +
  coord_flip() +
  theme_bw()
# ggsave('~/Dropbox/Frag Database (new)/analysis_apr19/figures/Sstd_forests_only.pdf', width = 290, height = 200, units = 'mm')
