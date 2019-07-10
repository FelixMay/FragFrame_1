# code to do prior predictive checks for fitting Bayesian models to fragment diversity data

# load the data
frag <- read_csv(paste0(path2data, '2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

# check lm fit
with(frag, lm(log(S_std_2) ~ c.lfs))

# get the fragment sizes for predicting slopes based on priors
x <- frag %>% 
  distinct(c.lfs) %>% 
  arrange(-desc(c.lfs))  

# visual inspection of priors for intercept (alpha) and slope (beta)
# recall we will model log(S) ~ c.lfs, i.e., both S and fragment size are log-transformed,
# and fragment size is mean centred too. This means that we want the predicted value to be
# ~mean(log(S_std_2)) when our predictor or x==0 (i.e., is the intercept of our linear model)
# And, we don't want to regularly predict impossibly strong (steep) relationships
alpha <- rnorm(50, 2.5, 0.2) # intercept ~ N(2.5, 0.2)
beta <- rnorm(50, 0, 0.5)  # slope ~ N(0, 1)

# simulate from these priors
y_pred <- tibble()
for(i in 1:length(beta)){
  temp = tibble(
    sim = i,
    x = x$c.lfs,
    y_temp = exp(alpha[i] + beta[i] * x))
  # join 'em together
  y_pred = bind_rows(y_pred, temp)
}


# slope inspection for predicting S_std_2: some a bit steep but not too bad
ggplot() + 
  geom_line(data = y_pred,
            aes(x = x, y = y_temp, group = sim),
            size = 0.1) +
  scale_y_continuous(trans = 'log',
                     name = 'y_pred') +
  geom_vline(xintercept = 0, lty =2) +
  geom_hline(data = frag %>% 
               summarise(mu = mean(S_std_2, na.rm=T)),
             aes(yintercept = mu), lty = 2) +
  geom_hline(data = frag %>% 
               summarise(max = max(S_std_2, na.rm=T)),
             aes(yintercept = max), lty = 2) +
  geom_hline(data = frag %>% 
               summarise(min = min(S_std_2, na.rm=T)),
             aes(yintercept = min), lty = 2) 


# what about N_std with these priors? No, intercept too low
ggplot() + 
  geom_line(data = y_pred,
            aes(x = x, y = y_temp, group = sim),
            size = 0.1) +
  scale_y_continuous(trans = 'log',
                     name = 'y_pred') +
  geom_vline(xintercept = 0, lty =2) +
  geom_hline(data = frag %>% 
               summarise(mu = mean(N_std, na.rm=T)),
             aes(yintercept = mu), lty = 2)  

with(frag, lm(log(N_std) ~ c.lfs))

alpha_N <- rnorm(100, 4, 0.2)
beta_N <- rnorm(100, 0, 0.5)

N_pred <- tibble()
for(i in 1:length(beta_N)){
  temp = tibble(
    sim = i,
    x = x$c.lfs,
    y_temp = exp(alpha_N[i] + beta_N[i] * x))
  # join 'em together
  N_pred = bind_rows(N_pred, temp)
}

# intercept a bit low still and a few strong slopes, but this looks ok
ggplot() + 
  geom_line(data = N_pred,
            aes(x = x, y = y_temp, group = sim),
            size = 0.1) +
  scale_y_continuous(trans = 'log',
                     name = 'N_pred') +
  geom_vline(xintercept = 0, lty =2) +
  geom_hline(data = frag %>% 
               summarise(mu = mean(N_std, na.rm=T)),
             aes(yintercept = mu), lty = 2) +
  geom_hline(data = frag %>% 
               summarise(mu = mean(N_std, na.rm=T)),
             aes(yintercept = mu), lty = 2) +
  geom_hline(data = frag %>% 
               summarise(max = max(N_std, na.rm=T)),
             aes(yintercept = max), lty = 2) +
  geom_hline(data = frag %>% 
               summarise(min = min(N_std, na.rm=T)),
             aes(yintercept = min), lty = 2)
