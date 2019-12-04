# source code from iNEXT package
#
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

# Set infinity values to NA
Inf_to_NA <- function(x)
{
   x <- as.numeric(x)
   x[!is.finite(x)] <- NA
   return(x)
}

# Get biodiversity indices from abundance vector
biodiv_abund <- function(abund_vec, N_std, n_base, cov_base){
   
   abund_vec <- as.numeric(abund_vec[abund_vec > 0])
   
   # get mobr indices
   mob <- calc_biodiv(abund_mat = matrix(abund_vec, nrow = 1, ncol = length(abund_vec)),
                      groups = "sample1",
                      index = c("S_n", "S_asymp", "S_PIE"), 
                      effort = n_base[1],
                      extrapolate = T,
                      return_NA = F)
   
   S_std2 <- rarefaction(as.numeric(abund_vec),
                         method = "indiv",
                         effort = N_std,
                         extrapolate = T,
                         quiet_mode = T)
   
   out <- c(S_std2    = as.numeric(S_std2),
            S_n       = mob$value[mob$index  == "S_n"],
            S_PIE     = mob$value[mob$index == "S_PIE"],
            S_cov     = NA,
            S_chao    = mob$value[mob$index == "S_asymp"]
   )
   
   S_cov <- try(estimateD(abund_vec, datatype = "abundance",
                                         base = "coverage",
                                         level = cov_base,
                                         conf = NULL))
   if (!is.error(S_cov)){
      out["S_cov"] <- S_cov[["q = 0"]]
   }
   
  return(out)
   
}

biodiv_resample <- function(dat, N_std, n_base, cov_base){
   sample1 <- sample(dat$species,
                     #size = round(N_std),
                     size = sum(dat$abundance),
                     prob = dat$abundance,
                     replace = T)
   abund_sample <- table(sample1)
   biodiv_sample <- biodiv_abund(abund_sample, N_std, n_base, cov_base)
}

biodiv_subplot <- function(abund_dat, N_std, n_base, cov_base, n_samples = 10){
   
   spec_abund <- data.frame(species  = names(abund_dat)[abund_dat > 0],
                            abundance = abund_dat[abund_dat > 0])
   
   #if (sum(spec_abund$abundance) == N_std){
      # biodiv <- biodiv_abund(spec_abund$abundance, n_base = n_base, cov_base = cov_base)
   #} else {
   
   # Calculate S_std by resample using N_std
   S_std1 <- replicate(n_samples, {
      sample1 <- sample(spec_abund$species,
                        size = N_std,
                        prob = spec_abund$abundance,
                        replace = T)
      abund_sample <- table(sample1)
      sum(abund_sample > 0)
      })
   
   samples <- replicate(n_samples, biodiv_resample(spec_abund,
                                                      N_std = N_std,
                                                      n_base = n_base,
                                                      cov_base = cov_base))
   samples <- rbind(S_std1, samples)
   
   biodiv_mean <- rowMeans(samples)
   biodiv_sd   <- apply(samples, 1, sd)
      
   names(biodiv_mean) <- paste(names(biodiv_mean), "mean", sep = "_")
   names(biodiv_sd)   <- paste(names(biodiv_sd), "sd", sep = "_")
   
   #}
      
   biodiv <- c(biodiv_mean, biodiv_sd)
   return(biodiv)
}

# https://stackoverflow.com/questions/23474729/convert-object-of-class-dist-into-data-frame-in-r
dist_to_dataframe <- function(inDist) {
   if (class(inDist) != "dist") stop("wrong input type")
   A <- attr(inDist, "Size")
   B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
   if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
   if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
   data.frame(
      row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
      col = rep(B[-length(B)], (length(B)-1):1),
      value = as.vector(inDist),
      stringsAsFactors = F
   )
}

# calculate all beta-diversity partitioning indices from sites x species table
# input: table with columns frag_id, frag_size_num, spec_1 .. spec_n

get_beta_part <- function(comat){
   
   comat <- arrange(comat, frag_size_num)
   class(comat) <- "data.frame"
   rownames(comat) <- comat$frag_id
   
   # pars <- expand.grid(coef = c("S","J","BS","BJ"),
   #                     quant = c(FALSE,TRUE), stringsAsFactors = F)
   
   # only calculate Jaccard and Ruzicka (= quantitative Jaccard) from Baselga family
   pars <- expand.grid(coef = "BJ",
                       quant = c(FALSE,TRUE),
                       stringsAsFactors = F)
   
   frag_pairs <- list()
   study_level <- list()
   
   for (i in 1:nrow(pars)){
      beta_part <- beta.div.comp(select(comat, -(frag_id:frag_size_num)),
                                 coef = pars$coef[i],
                                 quant = pars$quant[i])
      beta_tab <- dist_to_dataframe(beta_part$repl)
      names(beta_tab)[1:3] <- c("frag_x", "frag_y", "repl")
      
      beta_tab$rich <- dist_to_dataframe(beta_part$rich)[[3]]
     
      beta_tab$coef <- pars$coef[i]
      beta_tab$quant <- pars$quant[i]
      beta_tab$method <- beta_part$Note
      
      frag_pairs[[i]] <- beta_tab
      
      study_level[[i]] <- data.frame(coef = pars$coef[i],
                                     quant = pars$quant[i],
                                     prop_repl = beta_part$part["Repl/BDtotal"],
                                     prop_rich = beta_part$part["RichDif/BDtotal"],
                                     method =  beta_part$Note,
                                     stringsAsFactors = F, row.names = NULL)
   }
   
   return(list(fragments = bind_rows(frag_pairs), study = bind_rows(study_level)))
}

# resample fragment with given N
resample_fragment <- function(dat1, N){
   
   dat1 <- dat1 %>%
      select(frag_id, frag_size_num, species, abundance)
   
   dat_sample <- dat1[sample(1:nrow(dat1),
                             size = N,
                             prob = dat1$abundance,
                             replace = T),]
   dat_sample2 <- dat_sample %>%
      group_by(frag_id, frag_size_num, species) %>%
      count() %>%
      arrange(desc(n))
   names(dat_sample2) <- names(dat1)
   
   return(dat_sample2)
}

# get beta-diversity partitioning from bootstrap sample
boot_betapart <- function(index, dat_w){
   
   # randomly select one subplot per fragment
   dat_sample_w <- dat_w %>%
      select(-dataset_label, -sample_id, -frag_size_char, -sample_design) %>%
      group_by(frag_id) %>%
      sample_n(size = 1) %>%
      ungroup()
   
   # Calculate N, rel_sampling eff and N_std
   dat_sample_w$N <- rowSums(select(dat_sample_w, -(frag_id:sample_eff)))
   
   dat_sample_w <- dat_sample_w %>%
      mutate(rel_sample_eff = sample_eff/min(sample_eff),
             N_std = N/rel_sample_eff) %>%
      select(frag_id, frag_size_num, sample_eff, rel_sample_eff, N, N_std, everything())
   
   # convert to long format
   dat_sample_l <- dat_sample_w %>% 
      gather(key = species, value = abundance, -(frag_id:N_std)) %>%
      filter(abundance > 0)
   
   # resample sub-plots to N_std
   dat_sample_l_resample <- dat_sample_l %>%
      split(.$frag_id) %>%
      map2_dfr(round(dat_sample_w$N_std), resample_fragment)
   
   # convert to wide data
   dat_resample_w <- dat_sample_l_resample %>% 
      spread(key = species, value = abundance, fill = 0) 
   
   beta_part <- get_beta_part(dat_resample_w)
   
   return(beta_part)
}

# Function to calculate biodiversity indices for every fragment
get_biodiv <- function(data_set, n_thres = 5, fac_cont = 10,
                       method_abund = c("as_is","round","ceiling","multiply"),
                       n_resamples = 30){
   
   method_abund <- match.arg(method_abund)
   
   print(data_set$dataset_label[1])
   
   # test for non-integer abundances
   if (any(data_set$abundance %% 1 > 0)){
      print(paste("Non-integer abundances in ", data_set$dataset_label[1], sep = ""))
      
      if (method_abund == "round")
         data_set$abundance <- round(data_set$abundance)
      if (method_abund == "ceiling")
         data_set$abundance <- ceiling(data_set$abundance)
      if (method_abund == "multiply"){
         min_abund <- min(data_set$abundance[data_set$abundance > 0])
         data_set$abundance <- data_set$abundance / min_abund
      }
   }
   
   # impute fragment size for continuous fragments
   frag_cont <- data_set$frag_size_char == "continuous" & is.na(data_set$frag_size_num)
   data_set$frag_size_num[frag_cont] <- fac_cont * max(data_set$frag_size_num, na.rm = T)
   if (any(frag_cont))
      print(paste("Imputed area of continuous habitat in ", data_set$dataset_label [1], sep = ""))
   
   # sum up abundances in subplots
   # create wide data set
   dat_long <- data_set %>%
      group_by_at(vars(-abundance)) %>%
      summarise(abundance = sum(abundance)) %>%
      ungroup() %>%
      arrange(frag_size_num, sample_id)
   
   dat_wide <- dat_long %>%
      spread(key = species, value = abundance, fill = 0) %>%
      arrange(frag_size_num, sample_id)
   
   # abundance data only
   abund_mat <- as.data.frame(select(dat_wide, -(dataset_label:sample_design)))
   
   # one row per sub-sample
   dat_biodiv <- select(dat_wide, dataset_label:sample_design)
   
   # get total abundance and standardized total abundance per plot
   dat_biodiv <- dat_biodiv %>%
      mutate(rel_sample_eff = sample_eff/min(sample_eff),
             N = rowSums(abund_mat),
             N_std = N/rel_sample_eff,
             S_obs = rowSums(abund_mat > 0))
   
   # The following calculations are based on the raw data. Reference levels for n and
   # coverage might change, when we consider re-sampled data!
   
   # Determine base sample size for the sample-size-based rarefaction and extrapolation
   # Chao et al. 2014. Box 1
   n_max <- max(dat_biodiv$N)
   n_min_r <- min(2 * dat_biodiv$N)
   
   dat_biodiv$n_base <- round(min(n_max, n_min_r))   
   # dat_biodiv$n_base <- round(max(n_max, n_min_r))   

   # Chao et al. 2014 suggests using the maximum here,
   # we decided to go for a more conservative measure
   
   # Determine base coverage for the coverage-based rarefaction and extrapolation
   # Chao et al. 2014. Box 1
   
   # sample coverage
   dat_biodiv$coverage <- apply(abund_mat, 1,
                                function(x) {Chat.Ind(x, sum(x))} )
   
   # extrapolated coverage with sample size * 2
   cov_extra <- apply(abund_mat, 1,
                      function(x) {Chat.Ind(x, 2*sum(x))})
   
   # base coverage 
   dat_biodiv$cov_base <- min(max(dat_biodiv$coverage, na.rm = T),
                              min(cov_extra, na.rm = T))
   
   # calculate biodiv indices
   biodiv_indices <- mapply(biodiv_subplot,
                            split(abund_mat, 1:nrow(abund_mat)), 
                            N_std = round(dat_biodiv$N_std),
                            MoreArgs = list(n_base = dat_biodiv$n_base[1],
                                            cov_base = dat_biodiv$cov_base[1],
                                            n_samples = n_resamples)
                            )
   dat_biodiv <- cbind(dat_biodiv, data.frame(t(biodiv_indices)))
   
   # set indices to NA when there are no individuals
   empty_plots <- dat_biodiv$N == 0
   dat_biodiv$S_PIE_mean[empty_plots] <- NA
   dat_biodiv$S_PIE_sd[empty_plots] <- NA

   # set indices to NA when they are Inf or NaN
   dat_biodiv[,9:ncol(dat_biodiv)] <- lapply(dat_biodiv[,9:ncol(dat_biodiv)], Inf_to_NA)
   
   # set coverage standardized richness equal to observed
   # when observed coverage and base coverage equal 1
   # coverage_eq_1 <- dat_biodiv$coverage == 1 & dat_biodiv$cov_base == 1
   # dat_biodiv$S_cov_mean[coverage_eq_1] <- dat_biodiv$S_obs[coverage_eq_1]
   # What to do here with S_cov_sd
   
   # set S_n to NA when n_base < 5
   dat_biodiv$S_n_mean[dat_biodiv$n_base < n_thres] <- NA
   dat_biodiv$S_n_sd[dat_biodiv$n_base < n_thres] <- NA
   
   
   # average across sub-samples
   dat_biodiv_avg <- dat_biodiv %>%
      select(-sample_id, -sample_eff, -rel_sample_eff) %>%
      group_by_at(vars(dataset_label:sample_design)) %>%
      summarise_all(mean) %>%
      ungroup() %>%
      arrange(frag_size_num)
   
   # add mean and sd of largest fragment to all rows
   largest_frag <- dat_biodiv_avg %>% 
      filter(frag_size_num == max(frag_size_num)) %>%
      select(S_std1_mean:S_chao_sd) %>%
      summarise_all(mean)
   names(largest_frag) <- paste("exp", names(largest_frag), sep = "_")
   
   dat_out <- cbind(dat_biodiv_avg, largest_frag)
   
   dat_out <- dat_out %>%
      mutate(z_S_std  = (S_std1_mean - exp_S_std1_mean)/exp_S_std1_sd,
             z_S_n    = (S_n_mean - exp_S_n_mean)/exp_S_n_sd,
             z_S_PIE  = (S_PIE_mean - exp_S_PIE_mean)/exp_S_PIE_sd,
             z_S_cov  = (S_cov_mean - exp_S_cov_mean)/exp_S_cov_sd,
             z_S_chao = (S_chao_mean - exp_S_chao_mean)/exp_S_chao_sd
             ) 

   # beta-diversity partitioning ---------------------------------------------
   
   # analyse sampling design
   dat_design <- dat_long %>%
      select(frag_id, frag_size_num, sample_id, sample_eff) %>%
      distinct() %>%
      group_by(frag_id, frag_size_num) %>% 
      summarise(n_sub_samples = n(),
                range_sample_eff = max(sample_eff) - min(sample_eff)) %>%
      ungroup() %>%
      arrange(frag_size_num)
   
   # Equal sampling effort and one sub-sample per fragment
   if (all(dat_design$n_sub_samples == 1) && all(dat_design$range_sample_eff == 0)){

      beta_div_comp <- get_beta_part(select(dat_wide, -dataset_label, -sample_id,
                                                      -frag_size_char, -sample_eff, - sample_design))

      beta_part_frag <- beta_div_comp$fragments
      beta_part_study <- beta_div_comp$study
   
   } else { # use bootstrap approach
      
      beta_part_samples <- map(1:n_resamples, boot_betapart,
                               dat_w = dat_wide)
      
      beta_part_frag <- beta_part_samples %>% map_dfr("fragments")
      beta_part_study <- beta_part_samples %>% map_dfr("study")
   }
    
   # average across re-samples and partitioning methods
   beta_part_frag <- beta_part_frag %>%
      group_by(frag_x, frag_y, coef, quant, method) %>%
      summarise(repl = mean (repl),
                rich = mean (rich)) %>% ungroup()

   beta_part_study <- beta_part_study %>%
      group_by(coef, quant, method) %>%
      summarise(prop_repl = mean(prop_repl, na.rm = F),
                prop_rich = mean(prop_rich, na.rm = F)) %>% ungroup()

   # Add fragment areas difference and log-ratio
   dat_frag <- data_set %>%
      select(frag_id, frag_size_char, frag_size_num) %>%
      distinct()

   beta_part_frag <- beta_part_frag %>%
      left_join(dat_frag, by = c("frag_x" = "frag_id")) %>%
      left_join(dat_frag, by = c("frag_y" = "frag_id")) %>%
      mutate(diff_area = frag_size_num.x - frag_size_num.y,
             log10_ratio_area = log10(frag_size_num.x/frag_size_num.y)
      )

   # adjust beta-diversity partitioning dataframes
   beta_part_frag$dataset_label <- data_set$dataset_label[1]
   beta_part_frag$sample_design <- data_set$sample_design[1]

   beta_part_study$dataset_label <- data_set$dataset_label[1]
   beta_part_study$sample_design <- data_set$sample_design[1]

   beta_part_frag <- beta_part_frag %>%
      select(dataset_label,
             sample_design,
             frag_x,
             frag_y,
             frag_size_char.x,
             frag_size_num.x,
             frag_size_char.y,
             frag_size_num.y,
             diff_area,
             log10_ratio_area,
             everything())

   beta_part_study <- beta_part_study %>%
      select(dataset_label,
             sample_design,
             everything())

   out_list <- list("biodiv_frag" = dat_out,
                    "betapart_frag" = beta_part_frag,
                    "betapart_study" = beta_part_study
                    )
   
   return(out_list)
}    

# Execution of script -----------------------------------------------------

# read long format data file
infile <- path2Dropbox %+% "files_datapaper/Long_format_database/fragSAD_predicts_ewers.csv"
dat_all <- read.csv(infile, stringsAsFactors = F)
dim(dat_all)
str(dat_all)

head(dat_all)

frag_label <- dat_all %>% select(dataset_label, frag_id) %>% distinct()
dim(frag_label) 
data_set <- dat_all %>% filter(dataset_label == "Andresen_2003")
test <- get_biodiv(data_set)

data_set %>% 
   select(frag_id, sample_id, sample_eff) %>%
   distinct()

# get number of studies with specific sampling design
dat_all %>% 
   select(dataset_label, sample_design) %>%
   distinct() %>% count(sample_design)

# dat_all <- filter(dat_all, dataset_label != "Andresen_2003" & dataset_label != "Bernard_2007" )

parset <- expand.grid(fac_cont = c(2,10,100),
                      method_abund = c("as_is","round","ceiling","multiply"),
                      stringsAsFactors = F)
parset <- parset[c(1,2,3,8,11),]
#parset <- parset[2,]

for (i in 1:nrow(parset)){
   out1 <- dat_all %>%
      split(.$dataset_label) %>%
      map(get_biodiv,
          fac_cont = parset$fac_cont[i],
          method_abund = parset$method_abund[i],
          n_resamples = 200)
   
   out_biodiv_frag <- out1 %>% map_dfr("biodiv_frag")
   out_betapart_frag <- out1 %>% map_dfr("betapart_frag")
   out_betapart_study <- out1 %>% map_dfr("betapart_study")
   
   # prepare output
   outfile_name <- i %+% "_biodiv_frag_fcont_" %+% parset$fac_cont[i] %+%
      "_mabund_"%+% parset$method_abund[i] %+% ".csv"
   path2outfile <- path2Dropbox %+% "files_datapaper/Analysis/" %+% outfile_name
   write_csv(out_biodiv_frag, path2outfile)
   
   outfile_name <- i %+% "_betapart_frag_fcont_" %+% parset$fac_cont[i] %+%
      "_mabund_"%+% parset$method_abund[i] %+% ".csv"
   path2outfile <- path2Dropbox %+% "files_datapaper/Analysis/" %+% outfile_name
   write_csv(out_betapart_frag, path2outfile)

   outfile_name <- i %+% "_betapart_study_fcont_" %+% parset$fac_cont[i] %+%
      "_mabund_"%+% parset$method_abund[i] %+% ".csv"
   path2outfile <- path2Dropbox %+% "files_datapaper/Analysis/" %+% outfile_name
   write_csv(out_betapart_study, path2outfile)
}
   


# Note Problem with:
# Andresen 2003 --> mixed sampling design
# Bernard 2007 --> mixed sampling design
# Cadotte 2002 a --> check
