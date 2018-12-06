# code to fit models for fragmentation synthesis
# so far: bayesian framework for approximately ML-like results (i.e.,
# with non-informative priors - against the science of Gelman)

# 02/12/2018: to do:
# refit these models with more appropriate reference levels (time and matrix are inappropriate)


# load the data
# frag <- read_csv('~/Dropbox/Habitat loss meta-analysis/analysis/diversity_metadata.csv')

frag <- read_csv(paste(path2temp, "diversity_metadata.csv", sep = ""))


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
frag2$matrix.category <- factor(frag2$matrix.category, levels = c('light filter', 'medium filter', 'harsh filter'))
frag2$time.since.fragmentation <- factor(frag2$time.since.fragmentation,
                                         levels = c("recent (less than 20 years)", "intermediate (20-100 years)", "long (100+ years)"))

# fit two-way interaction between fragment size and: taxa, veg.fragment, matrix category and time since fragmentation
# S_std
lS_std_fS_taxa <- brm(lSstd ~ -1 + c.lfs*taxa + (c.lfs | filename), 
                   data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                     filter(!is.na(taxa)),
                   cores = 4, chains = 4)

lS_std_fS_veg <- brm(lSstd ~ c.lfs*veg.fragment + (c.lfs | filename), 
                  data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                    # leave out the one mangrove study, and those without data
                    filter(veg.fragment!='mangrove' & !is.na(veg.fragment)),
                  cores = 4, chains = 4)

lS_std_fS_matrix <- brm(lSstd ~ c.lfs*matrix.category + (c.lfs | filename), 
                     data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                       # leave out studies without data
                       filter(!is.na(matrix.category)),
                     cores = 4, chains = 4)

lS_std_fS_time <- brm(lSstd ~ c.lfs*time.since.fragmentation + (c.lfs | filename), 
                   data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                     # leave out studies without data
                     filter(!is.na(time.since.fragmentation)),
                   cores = 4, chains = 4)

lS_std_fS_2way_all <- brm(lSstd ~ c.lfs + taxa + veg.fragment + matrix.category + time.since.fragmentation +
                            c.lfs:taxa + c.lfs:veg.fragment + c.lfs:matrix.category + c.lfs:time.since.fragmentation + 
                            (c.lfs | filename), 
                      data = frag2 %>% filter(!filename %in% problems$filename & veg.fragment!='mangrove') %>% 
                        # leave out studies without data
                        filter(!is.na(taxa) & !is.na(veg.fragment) & 
                                 !is.na(matrix.category) & !is.na(time.since.fragmentation)),
                      cores = 4, chains = 4)

# S_n
lSn_fS_taxa <- brm(lSn ~ c.lfs*taxa + (c.lfs | filename), 
                   data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                     filter(!is.na(taxa)),
                   cores = 4, chains = 4)

lSn_fS_veg <- brm(lSn ~ c.lfs*veg.fragment + (c.lfs | filename), 
                  data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                    # leave out the one mangrove study, and those without data
                    filter(veg.fragment!='mangrove' & !is.na(veg.fragment)),
                  cores = 4, chains = 4)

lSn_fS_matrix <- brm(lSn ~ c.lfs*matrix.category + (c.lfs | filename), 
                     data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                       # leave out studies without data
                       filter(!is.na(matrix.category)),
                     cores = 4, chains = 4)

lSn_fS_time <- brm(lSn ~ c.lfs*time.since.fragmentation + (c.lfs | filename), 
                   data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                     # leave out studies without data
                     filter(!is.na(time.since.fragmentation)),
                   cores = 4, chains = 4)

lSn_fS_2way_all <- brm(lSn ~ c.lfs + taxa + veg.fragment + matrix.category + time.since.fragmentation +
                            c.lfs:taxa + c.lfs:veg.fragment + c.lfs:matrix.category + c.lfs:time.since.fragmentation + 
                            (c.lfs | filename), 
                          data = frag2 %>% filter(!filename %in% problems$filename & veg.fragment!='mangrove') %>% 
                            # leave out studies without data
                            filter(!is.na(taxa) & !is.na(veg.fragment) & 
                                     !is.na(matrix.category) & !is.na(time.since.fragmentation)),
                          cores = 4, chains = 4)

# repeat for S_PIE
lS_PIE_fS_taxa <- brm(lS_PIE ~ c.lfs*taxa + (c.lfs | filename), 
                   data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                     # leave out studies without data
                     filter(!is.na(taxa)),
                   cores = 4, chains = 4)

lS_PIE_fS_veg <- brm(lS_PIE ~ c.lfs*veg.fragment + (c.lfs | filename), 
                  data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                    # leave out the one mangrove study, and those without data
                    filter(veg.fragment!='mangrove' & !is.na(veg.fragment)),
                  cores = 4, chains = 4)

lS_PIE_fS_matrix <- brm(lS_PIE ~ c.lfs*matrix.category + (c.lfs | filename), 
                     data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                       # leave out studies without data
                       filter(!is.na(matrix.category)),
                     cores = 4, chains = 4)

lS_PIE_fS_time <- brm(lS_PIE ~ c.lfs*time.since.fragmentation + (c.lfs | filename), 
                   data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                     # leave out studies without data
                     filter(!is.na(time.since.fragmentation)),
                   cores = 4, chains = 4)

lS_PIE_fS_2way_all <- brm(lS_PIE ~ c.lfs + taxa + veg.fragment + matrix.category + time.since.fragmentation +
                         c.lfs:taxa + c.lfs:veg.fragment + c.lfs:matrix.category + c.lfs:time.since.fragmentation + 
                         (c.lfs | filename), 
                       data = frag2 %>% filter(!filename %in% problems$filename & veg.fragment!='mangrove') %>% 
                         # leave out studies without data
                         filter(!is.na(taxa) & !is.na(veg.fragment) & 
                                  !is.na(matrix.category) & !is.na(time.since.fragmentation)),
                       cores = 4, chains = 4)

# repeat for N_std
Nstd_lognorm_fS_taxa <- brm(N_std ~ c.lfs*taxa + (c.lfs | filename), 
                      data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                        # leave out studies without data
                        filter(!is.na(taxa) & N_std > 0),
                      family = 'lognormal',
                      cores = 4, chains = 4)

Nstd_lognorm_fS_veg <- brm(N_std ~ c.lfs*veg.fragment + (c.lfs | filename), 
                     data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                       # leave out the one mangrove study, and those without data
                       filter(veg.fragment!='mangrove' & !is.na(veg.fragment) & N_std > 0),
                     family = 'lognormal',
                     cores = 4, chains = 4)

Nstd_lognorm_fS_matrix <- brm(N_std ~ c.lfs*matrix.category + (c.lfs | filename), 
                        data = frag2 %>% filter(!filename %in% problems$filename & N_std > 0) %>% 
                          # leave out studies without data
                          filter(!is.na(matrix.category)),
                        family = 'lognormal',
                        cores = 4, chains = 4)

Nstd_lognorm_fS_time <- brm(N_std ~ c.lfs*time.since.fragmentation + (c.lfs | filename), 
                      data = frag2 %>% filter(!filename %in% problems$filename & N_std > 0) %>% 
                        # leave out studies without data
                        filter(!is.na(time.since.fragmentation)),
                      family = 'lognormal',
                      cores = 4, chains = 4)

Nstd_lognorm_fS_2way_all <- brm(N_std ~ c.lfs + taxa + veg.fragment + matrix.category + time.since.fragmentation +
                            c.lfs:taxa + c.lfs:veg.fragment + c.lfs:matrix.category + c.lfs:time.since.fragmentation + 
                            (c.lfs | filename), 
                          data = frag2 %>% filter(!filename %in% problems$filename & veg.fragment!='mangrove' & N_std > 0) %>% 
                            # leave out studies without data
                            filter(!is.na(taxa) & !is.na(veg.fragment) & 
                                     !is.na(matrix.category) & !is.na(time.since.fragmentation)),
                          family = 'lognormal',
                          cores = 4, chains = 4)

#refit N_std models with log-transformed response for compatibility with other models
# repeat for N_std
lNstd_fS_taxa <- brm(lNstd ~ c.lfs*taxa + (c.lfs | filename), 
                            data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                              # leave out studies without data
                              filter(!is.na(taxa) & N_std > 0),
                            cores = 4, chains = 4)

lNstd_fS_veg <- brm(lNstd ~ c.lfs*veg.fragment + (c.lfs | filename), 
                           data = frag2 %>% filter(!filename %in% problems$filename) %>% 
                             # leave out the one mangrove study, and those without data
                             filter(veg.fragment!='mangrove' & !is.na(veg.fragment) & N_std > 0),
                           cores = 4, chains = 4)

lNstd_fS_matrix <- brm(lNstd ~ c.lfs*matrix.category + (c.lfs | filename), 
                              data = frag2 %>% filter(!filename %in% problems$filename & N_std > 0) %>% 
                                # leave out studies without data
                                filter(!is.na(matrix.category)),
                              cores = 4, chains = 4)

lNstd_fS_time <- brm(lNstd ~ c.lfs*time.since.fragmentation + (c.lfs | filename), 
                            data = frag2 %>% filter(!filename %in% problems$filename & N_std > 0) %>% 
                              # leave out studies without data
                              filter(!is.na(time.since.fragmentation)),
                            cores = 4, chains = 4)

setwd('~/Dropbox/1current/fragmentation_synthesis/results/')
save(lS_std_fS_taxa, lS_std_fS_veg, lS_std_fS_matrix, lS_std_fS_time, lSn_fS_2way_all,
     lSn_fS_taxa, lSn_fS_veg, lSn_fS_matrix, lSn_fS_time, lSn_fS_2way_all,
     lS_PIE_fS_taxa, lS_PIE_fS_veg, lS_PIE_fS_matrix, lS_PIE_fS_time, lS_PIE_fS_2way_all,
     Nstd_lognorm_fS_taxa, Nstd_lognorm_fS_veg, Nstd_lognorm_fS_matrix, Nstd_lognorm_fS_time, Nstd_lognorm_fS_2way_all,
     lNstd_fS_taxa, lNstd_fS_veg, lNstd_fS_matrix, lNstd_fS_time,
     file = 'two_interactions_brms_fits.Rdata')

##----load model fits-----
load('~/Dropbox/1current/fragmentation_synthesis/results/two_interactions_brms_fits.Rdata')

# coefficient plots
two_way_fixed_coefs <- bind_rows(fixef(lS_std_fS_taxa) %>% 
                                as_tibble() %>% 
                                mutate(response = 'S_std',
                                       model = 'fS_x_taxa',
                                       coef = rownames(fixef(lS_std_fS_taxa))),
                              fixef(lS_std_fS_matrix) %>% 
                                as_tibble() %>% 
                                mutate(response = 'S_std',
                                       model = 'fS_x_matrix',
                                       coef = rownames(fixef(lS_std_fS_matrix))),
                              fixef(lS_std_fS_time) %>% 
                                as_tibble() %>% 
                                mutate(response = 'S_std',
                                       model = 'fS_x_time',
                                       coef = rownames(fixef(lS_std_fS_time))),
                              fixef(lS_std_fS_veg) %>% 
                                as_tibble() %>% 
                                mutate(response = 'S_std',
                                       model = 'fS_x_vegFrag',
                                       coef = rownames(fixef(lS_std_fS_veg))),
                              # Sn
                              fixef(lSn_fS_taxa) %>% 
                                as_tibble() %>% 
                                mutate(response = 'Sn',
                                       model = 'fS_x_taxa',
                                       coef = rownames(fixef(lSn_fS_taxa))),
                              fixef(lSn_fS_matrix) %>% 
                                as_tibble() %>% 
                                mutate(response = 'Sn',
                                       model = 'fS_x_matrix',
                                       coef = rownames(fixef(lSn_fS_matrix))),
                              fixef(lSn_fS_time) %>% 
                                as_tibble() %>% 
                                mutate(response = 'Sn',
                                       model = 'fS_x_time',
                                       coef = rownames(fixef(lSn_fS_time))),
                              fixef(lSn_fS_veg) %>% 
                                as_tibble() %>% 
                                mutate(response = 'Sn',
                                       model = 'fS_x_vegFrag',
                                       coef = rownames(fixef(lSn_fS_veg))),
                              # S_PIE
                              fixef(lS_PIE_fS_taxa) %>% 
                                as_tibble() %>% 
                                mutate(response = 'S_PIE',
                                       model = 'fS_x_taxa',
                                       coef = rownames(fixef(lS_PIE_fS_taxa))),
                              fixef(lS_PIE_fS_matrix) %>% 
                                as_tibble() %>% 
                                mutate(response = 'S_PIE',
                                       model = 'fS_x_matrix',
                                       coef = rownames(fixef(lS_PIE_fS_matrix))),
                              fixef(lS_PIE_fS_time) %>% 
                                as_tibble() %>% 
                                mutate(response = 'S_PIE',
                                       model = 'fS_x_time',
                                       coef = rownames(fixef(lS_PIE_fS_time))),
                              fixef(lS_PIE_fS_veg) %>% 
                                as_tibble() %>% 
                                mutate(response = 'S_PIE',
                                       model = 'fS_x_vegFrag',
                                       coef = rownames(fixef(lS_PIE_fS_veg))),
                              # N_std
                              fixef(lNstd_fS_taxa) %>% 
                                as_tibble() %>% 
                                mutate(response = 'N_std',
                                       model = 'fS_x_taxa',
                                       coef = rownames(fixef(Nstd_lognorm_fS_taxa))),
                              fixef(lNstd_fS_matrix) %>% 
                                as_tibble() %>% 
                                mutate(response = 'N_std',
                                       model = 'fS_x_matrix',
                                       coef = rownames(fixef(Nstd_lognorm_fS_matrix))),
                              fixef(lNstd_fS_time) %>% 
                                as_tibble() %>% 
                                mutate(response = 'N_std',
                                       model = 'fS_x_time',
                                       coef = rownames(fixef(Nstd_lognorm_fS_time))),
                              fixef(lNstd_fS_veg) %>% 
                                as_tibble() %>% 
                                mutate(response = 'N_std',
                                       model = 'fS_x_vegFrag',
                                       coef = rownames(fixef(Nstd_lognorm_fS_veg)))) %>% 
  # add indicator for CIs that differ from zero
  mutate(overlap_0 = ifelse((Q2.5 < 0 & Q97.5 <0) |
                              (Q2.5 > 0 & Q97.5 > 0),
                            '1', '0'),
         labels = coef)
# fix labels for plotting
two_way_fixed_coefs$model <- factor(two_way_fixed_coefs$model,
                                    levels = c("fS_x_taxa", "fS_x_matrix", "fS_x_time", "fS_x_vegFrag"),
                                    labels = c('Frag. size x taxa', 'Frag. size x matrix harshness', 'Frag. size x time since frag.', 'Frag. size x frag. vegetation'))
two_way_fixed_coefs$labels <- factor(two_way_fixed_coefs$labels,
                                     levels = c('Intercept', 'c.lfs',
                                                'taxabirds', 'taxainvertebrates', 'taxamammals', 'taxaplants',
                                                'veg.fragmentgrassland', 'veg.fragmentshrublandDsteppe',
                                                'time.since.fragmentationintermediate20M100years',  'time.since.fragmentationlong100Pyears', 
                                                'c.lfs:taxabirds', 'c.lfs:taxainvertebrates', 'c.lfs:taxamammals', 'c.lfs:taxaplants', 
                                                'matrix.categorymediumfilter', 'matrix.categoryharshfilter',
                                                'c.lfs:matrix.categorymediumfilter',  'c.lfs:matrix.categoryharshfilter', 
                                                'c.lfs:time.since.fragmentationintermediate20M100years', 'c.lfs:time.since.fragmentationlong100Pyears',
                                                'c.lfs:veg.fragmentgrassland', 'c.lfs:veg.fragmentshrublandDsteppe'),
                                     labels = c('Intercept', 'Fragment size (log)',
                                                'Birds', 'Invertebrates', 'Mammals', 'Plants',
                                                'Grassland', 'Shrubland / steppe',
                                                'Time since fragmentation (20-100 years)',  'Time since fragmentation (>100 years)', 
                                                'Frag. size x birds', 'Frag. size x invertebrates', 'Frag. size x mammals', 'Frag. size x plants', 
                                                'Medium matrix', 'Harsh matrix',
                                                'Frag. size x medium matrix',  'Frag. size x harsh matrix', 
                                                'Frag. size x Time since fragmentation (20-100 years)',  'Frag. size x Time since fragmentation (>100 years)', 
                                                'Frag size x grassland', 'Frag. size x Shrubland / steppe'))


setwd('~/Dropbox/Habitat loss meta-analysis/analysis/figs/')
ggplot() +
  facet_wrap(~model, scales = 'free') +
  geom_point(data = two_way_fixed_coefs,
             aes(x = labels, y = Estimate, group = response,
                 colour = response, shape = overlap_0),
             position = position_dodge(width = 0.8),
             size = 2) +
  geom_errorbar(data = two_way_fixed_coefs,
                aes(x = labels, ymin = Q2.5, ymax = Q97.5,
                    colour = response, group = response, linetype = overlap_0),
                width = 0,
                position = position_dodge(width = 0.8)) +
  scale_shape_manual(guide = F, values = c('1' = 19, '0' = 1)) +
  scale_linetype_manual(guide = F, values = c('1' = 1, '0' = 2)) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = 'Coefficient', 
       y = 'Coefficient estimate',
       caption = 'Coefficient estimates for models with two-way interactions. 
       Intercept and Frag. size effect in all panels is for the reference level of the term:
       taxa = amphibians/reptiles, matrix harshness = light,
       time since frag. = short (< 20 years), frag. vegetation = forest.
       Filled symbols and solid lines represent coefficients that differ from zero.') +
  coord_flip() +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0))

# ggsave('two_way_interactions.pdf', width = 330, height = 200, units = 'mm')
