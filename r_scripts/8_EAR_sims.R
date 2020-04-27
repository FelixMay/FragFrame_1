# code to simulate species loss from EAR with random sampling and ecosystem decay

# the total area (before fragmentation)
A = 1e4
# total species richness in A
S = 1e3
# fragment sizes of smaller areas than A
a = seq(from = 1, to = A, length.out = 200)
#  even SAD
cv_abund=1
# total community size (total number of individuals of all species combined)
totalN = 1e4

# generate regional community; log-normal SAD, possion 
comm <- sim_poisson_community(S, totalN, 'lnorm',
                              xrange = c(0,sqrt(A)),
                              yrange = c(0,sqrt(A)))
# plot(comm)

# we need n_i, the number of individuals for each species
temp <- table(comm$census$species)
df <- tibble(species = names(temp),
             n_i = temp)

# recall a = subarea lost, A = total area

# Estimate species remaining as random subsample from A using S_ear = (a/A)^n_i, summed from i to S species
S_ear <- matrix(data = NA,
                nrow = length(a),
                ncol = 3)

# define beta (here, equal to our empirical slope estimate) for decay
beta = 0.06

for(i in 1:length(a)){
  prop_lost =  a[i]/A
  area_remaining = A - a[i]
  n_i = df$n_i
  # calculate S_random = S_ear (i.e., # species in A-a)
  S_random = sum(prop_lost ^ n_i)
  # This is where I've changed things from the initial submission:
  # we expected to lose additional species a under ecosystem decay, 
  # I haven't derived the logic as we did for the initial submission,
  # but the following, where I multiply our slope estimate by the size of the area lost gives a pattern
  # similar to that in Felix's new conceptual diagram
  S_decay = S_random + beta*a[i]
  S_ear[i,] = c(a[i], S_random, S_decay)
}

colnames(S_ear) <- c('a', 'S_random', 'S_decay')


# to tibble for plotting
S_ear <- as_tibble(S_ear) %>% 
  mutate(S = S,
         A = A)

# plot as number of extinctions as a function of area lost (as per Felix's new figure)
ggplot() +
  geom_line(data = S_ear,
            aes(x = a, y = S_random, colour = 'S_random'),
            size = 2) +
  geom_line(data = S_ear %>% filter(S_decay < 1001),
            aes(x = a, y = S_decay, colour = 'S_decay'),
            size = 2) +
  scale_color_manual(name = '',
                     values = c('S_random' = '#7570b3',
                                'S_decay' = '#d95f02'),
                     labels = c(expression(S[decay]), expression(S[passive]))) +
  labs(x = 'Habitat area lost',
       y = 'Number of extinctions') +
  theme_classic() +
  theme(legend.position = c(0.1,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        text = element_text(size = 7))

# 1.5 column width
ggsave(paste0(path2wd, 'extended_data_gis_tabs/Ex_Dat_Fig9.png'),
       width = 120, height = 120, units = 'mm')