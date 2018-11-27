# code to fit preliminary models for fragmentation synthesis
# preliminary models will use ML, but I plan to use Bayesian methods for 
# final analysis
library(lme4)

# load the data
frag <- read_csv('~/Dropbox/Habitat loss meta-analysis/analysis/diversity_metadata.csv')

# remove observations without fragment size
frag2 <- frag %>% 
  filter(!is.na(entity.size))

##--create some covariates for easier workflow--
# mean-centred log(fragment.size)
frag2$c.lfs <- log(frag2$entity.size) - mean(log(frag2$entity.size))
frag2$lSn <- log(frag2$S_n)

# simplest model: richness as a function of fragment size; allow fragment size to vary by study
summary(m1_fragSize <- lmer(lSn ~ c.lfs + (c.lfs | filename), 
                             data = frag2))

# inspect model fit
m1_inspect <- broom::augment(m1_fragSize) %>% 
  as_tibble() %>% 
  # add some other covariates that were not modelled here
  inner_join(frag2 %>% distinct(filename, time.since.fragmentation,
                              taxa, veg.fragment, matrix.category, ratio.min.max.fragment.size2),
             by = 'filename' )



ggplot() +
  geom_boxplot(data = m1_inspect,
               aes(filename, .resid)) +
  coord_flip()

ggplot() +
  geom_boxplot(data = m1_inspect,
               aes(veg.fragment, .resid)) +
  coord_flip()
ggplot() +
  geom_boxplot(data = m1_inspect,
               aes(matrix.category, .resid)) +
  coord_flip()

ggplot() +
  geom_point(data = m1_inspect,
               aes(c.lfs, .resid))

ggplot() +
  geom_point(data = m1_inspect,
             aes(time.since.fragmentation, .resid))

ggplot() +
  geom_point(data = m1_inspect,
             aes(ratio.min.max.fragment.size2, .resid))

ggplot() +
  geom_boxplot(data = m1_inspect %>% filter(!is.na(taxa)),
    aes(taxa, .resid))
        

summary(m2_fS_taxa <- lmer(lSn ~ -1 + c.lfs*taxa + (c.lfs | filename), 
                            data = frag2, REML = F))
car::Anova(m2_fS_taxa)
# alternate for taxa
summary(m2_fS_taxa_2 <- lmer(lSn ~ -1 + c.lfs*taxa + (c.lfs | taxa/filename),
                             data = frag2, REML = F))
AIC(m2_fS_taxa, m2_fS_taxa_2) # no support for adding taxa as grouping covariate

summary(m3_fS_taxa_matrix <- lmer(lSn ~ -1 + c.lfs*taxa*matrix.category + (c.lfs | filename), 
                           data = frag2, REML = T))
car::Anova(m3_fS_taxa_matrix)

summary(m4_fS_taxa_matrix <- lmer(lSn ~ -1 + c.lfs*taxa*matrix.category*veg.fragment + (c.lfs | filename), 
                                  data = frag2, REML = T))
car::Anova(m4_fS_taxa_matrix) # looks like there might be support for a fragmentSize x matrixCategory x vegFragment interaction
                              # and, fragmentSize x taxa,

# most complex
summary(m5_fS_taxa_matrix <- lmer(lSn ~ -1 + c.lfs*taxa*matrix.category*veg.fragment*time.since.fragmentation + (c.lfs | filename), 
                                  data = frag2, REML = T))
car::Anova(m5_fS_taxa_matrix) # looks like there might be support for a fragmentSize x matrixCategory x vegFragment interaction
# and, fragmentSize x taxa,