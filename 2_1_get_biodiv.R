# source code from iNEXT package
###############################################
# Abundance-based sample coverage
# 
# \code{Chat.Ind} Estimation of abundance-based sample coverage function
# 
# @param x a vector of species abundances
# @param m a integer vector of rarefaction/extrapolation sample size
# @return a vector of estimated sample coverage function
# @export
Chat.Ind <- function(x, m){
   x <- x[x>0]
   n <- sum(x)
   f1 <- sum(x == 1)
   f2 <- sum(x == 2)
   f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
   A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
   Sub <- function(m){
      #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
      if(m < n) {
         xx <- x[(n-x)>=m]
         out <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
      }
      if(m == n) out <- 1-f1/n*A
      if(m > n) out <- 1-f1/n*A^(m-n+1)
      out
   }
   sapply(m, Sub)		
}

######
Inf_to_NA <- function(x)
{
   x <- as.numeric(x)
   x[!is.finite(x)] <- NA
   return(x)
}




############################################
# Function to calculate biodiversity indices for every fragment

get_biodiv <- function(data_set), n_thres = 5){
   
   # sum abundances in same fragments
   dat_abund <- data_set %>%
      group_by(frag_id, species) %>%
      summarise(abundance = sum(abundance))
   
   dat_wide <- dat_abund %>% spread(key = species, value = abundance, fill = 0)
   
   # prepare output data
   dat_sample <- data_set %>%
      select(-species, -abundance) %>%
      distinct()
   
   dat_sample_eff <- dat_sample %>%
      group_by(frag_id) %>%
      summarise(sample_eff = sum(sample_eff))
   
   dat_frag <- dat_sample %>% 
      select(-sample_id, -sample_eff) %>%
      distinct() %>%
      left_join(dat_sample_eff)
   
   # get biodiversity indices
   dat_biodiv <- data.frame(frag_id = dat_wide$frag_id,
                            N = rowSums(dat_wide[,-1]))
   
   # Determine base sample size for the sample-size-based rarefaction and extrapolation
   # Chao et al. 2014. Box 1
   n_max <- max(dat_biodiv$N)
   n_min_r <- min(2 * dat_biodiv$N)
   
   dat_biodiv$n_base <- min(n_max, n_min_r)   

   # Chao et al. 2014 suggests using the maximum here,
   # we decided to go for a more conservative measure
   
   # Determine base coverage for the coverage-based rarefaction and extrapolation
   # Chao et al. 2014. Box 1
   
   # sample coverage
   dat_biodiv$coverage <- apply(dat_wide[,-1], 1,
                                function(x) {Chat.Ind(x, sum(x))} )
   
   # extrapolated coverage with sample size * 2
   cov_extra <- apply(dat_biodiv[,-1], 1,
                      function(x) {Chat.Ind(x, 2*sum(x))})
   
   # base coverage 
   dat_biodiv$cov_base <- min(max(dat_biodiv$coverage, na.rm = T),
                              min(cov_extra, na.rm = T))
    
   # get mobr indices
   mob <- calc_biodiv(abund_mat = as.data.frame(dat_wide[,-1]),
                      groups = rep("frag", nrow(dat_wide)),
                      index = c("S", "S_n", "S_PIE"), 
                      effort = dat_biodiv$n_base[1],
                      extrapolate = T,
                      return_NA = F)
   
   
   # observed richness
   dat_biodiv$S_obs <- mob$value[mob$index == "S"]
   
   # rarefied richness
   if (dat_biodiv$n_base[1] >= n_thres){
      dat_biodiv$S_n <- mob$value[mob$index  == "S_n"]
   } else {
      dat_biodiv$S_n <- NA
   }
   
   # ENS PIE
   dat_biodiv$S_PIE <- mob$value[mob$index == "S_PIE"]
   
   # coverage standardized richness
   div_indi$S_cov <- rep(NA, nrow(div_indi)) 
   
   S_cov_std <- lapply(data.frame(t(dat_wide[,-1])),
                       function(x) try(estimateD(x, datatype = "abundance",
                                                 base = "coverage",
                                                 level = dat_biodiv$cov_base[1],
                                                 conf = NULL)))
   succeeded <- !sapply(S_cov_std, is.error)
   
   if (sum(succeeded) > 0){
      dat_biodiv$S_cov[succeeded] <- sapply(S_cov_std[succeeded],"[[","q = 0")
   }
   
   
   #######################################################################
   
   ### simple diversity indices
   
   # standardize sampling effort
   # smallest sampling unit:  sampling effort = 1
   rel_sample_eff <- div_indi$sample_effort/min(div_indi$sample_effort)
   
 
   
   # n_min <- min(div_indi$N)
   
   # No. of individuals standardized by sampling effort
   
   #div_indi$N_std <-  div_indi$N/div_indi$sample_effort ## standardize abundance by absolute sampling effort
   div_indi$N_std <-  div_indi$N/rel_sample_eff ## standardize abundance by relative sampling effort
   
   # richness standardized by mean sampling effort in smallest sampling unit
   div_indi$S_std <- NA
   for (i in 1:nrow(div_indi)){
      if (div_indi$N_std[i] >= n_thres)   
         div_indi$S_std[i] <- rarefaction(dat_abund_pool2[,i],
                                          method = "indiv",
                                          effort = round(div_indi$N_std[i]))        
   }

  
   
   # div_indi$ENS_pie <- diversity(t(dat_abund_pool2), index = "invsimpson")
   
   #####################################################################
   # I do not remember when and how we defined the following calulation
   # S_plot <- colSums(dat_abund > 0)    # observed species richness per plot before plots within the same fragments have been pooled
   # div_indi$S_plot <- ifelse(div_indi$sample_design == "pooled", NA, aggregate(S_plot, 
   #                          by = list(dat_head_t$entity.id,
   #                                    dat_head_t$entity.size.rank),
   #                          FUN = mean)$x) # mean of plot level species richness, comparison only meaningful 
   # if sampling effort among plots is equal, i.e. sampling_design either standardized or
   # multiple standardized subsamples within fragments
   #####################################################################
   
  
   
   # I am not sure these three lines are neede?
   # cov_eq_1 <- div_indi$coverage >= 1.0
   # cov_eq_1[is.na(cov_eq_1)] <- FALSE
   # div_indi$D0_hat[cov_eq_1] <- div_indi$S_obs[cov_eq_1]

   # set indices to NA when there are no individuals
   empty_plots <- div_indi$N == 0
   div_indi$S_asymp[empty_plots] <- NA
   div_indi$S_PIE[empty_plots] <- NA

   
   # set indices to NA when they are Inf or NaN
   div_indi[,9:ncol(div_indi)] <- lapply(div_indi[,9:ncol(div_indi)], Inf_to_NA)
   
   return(div_indi)
}    

################################################################################
# execution of script

# read long format data file

infile <- path2Dropbox %+% "files_datapaper/Long_format_database/fragSAD_long.csv"
dat_long <- read.csv(infile, stringsAsFactors = F)
dim(dat_long)
str(dat_long)

head(dat_long)

data_set <- dat_long %>% filter(dataset_id == "53")



