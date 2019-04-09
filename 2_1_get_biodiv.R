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
get_beta_part <- function(comat){
   
   comat <- arrange(comat, frag_size_num)
   class(comat) <- "data.frame"
   rownames(comat) <- comat$frag_id
   
   pars <- expand.grid(coef = c("S","J","BS","BJ"),
                       quant = c(FALSE,TRUE), stringsAsFactors = F)
   
   frag_pairs <- list()
   study_level <- list()
   
   for (i in 1:nrow(pars)){
      beta_part <- beta.div.comp(comat[,-(1:3)],
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

# Sample subplots from community matrix
betapart_subplots <- function(index, dat1, n){
   dat_sample <- dat1 %>%
      group_by(frag_id) %>%
      sample_n(n)
   
   dat_frag <- dat_sample %>%
      select(frag_id, frag_size_num) %>%
      distinct
   
   # sum abundances in the same fragment
   dat_sample2 <- dat_sample %>%
      select(-frag_size_num, -sample_id) %>%
      group_by(frag_id) %>%
      summarise_all(sum)
   
   dat_sample3 <- dat_frag %>%
      left_join(dat_sample2) %>%
      arrange(frag_size_num)
   
   beta_part <- get_beta_part(dat_sample3)
   return(beta_part)
}

# Resample abundance data and calculate beta-div partitioning
betapart_resample <- function(index, dat1, N_std){
   
   dat_abund_resample <- dat1 %>%
      split(.$frag_id) %>%
      map2_dfr(N_std, resample_fragment)
   
   dat_wide_resample <- dat_abund_resample %>% 
      spread(key = species, value = abundance, fill = 0) 
   
   beta_part <- get_beta_part(dat_wide_resample)
   return(beta_part)
}

# resample fragment with given N
resample_fragment <- function(dat1, N){
   dat_sample <- dat1[sample(1:nrow(dat1),
                             size = N,
                             prob = dat1$abundance,
                             replace = T), 1:3]
   dat_sample2 <- dat_sample %>%
      group_by(frag_id, frag_size_num, species) %>%
      count() %>%
      arrange(desc(n))
   names(dat_sample2) <- names(dat1)
   
   return(dat_sample2)
}


# Function to calculate biodiversity indices for every fragment
get_biodiv <- function(data_set, n_thres = 5, fac_cont = 10,
                       method_abund = c("as_is","round","ceiling","multiply"),
                       n_resamples = 100){
   
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
   
   # sum abundances in same fragments
   dat_abund <- data_set %>%
      group_by(frag_id, frag_size_num, species) %>%
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
      left_join(dat_sample_eff) %>%
      arrange(frag_size_num)
   
   dat_wide <- dat_frag %>%
      select(frag_id, sample_eff) %>%
      left_join(dat_wide) %>%
      arrange(frag_size_num)
   
   # get biodiversity indices ----
   dat_biodiv <- data.frame(frag_id = dat_wide$frag_id,
                            sample_eff = dat_wide$sample_eff,
                            N = rowSums(dat_wide[,-(1:3)]),
                            stringsAsFactors = F)
   
   # Determine base sample size for the sample-size-based rarefaction and extrapolation
   # Chao et al. 2014. Box 1
   n_max <- max(dat_biodiv$N)
   n_min_r <- min(2 * dat_biodiv$N)
   
   dat_biodiv$n_base <- round(min(n_max, n_min_r))   

   # Chao et al. 2014 suggests using the maximum here,
   # we decided to go for a more conservative measure
   
   # Determine base coverage for the coverage-based rarefaction and extrapolation
   # Chao et al. 2014. Box 1
   
   # sample coverage
   dat_biodiv$coverage <- apply(dat_wide[,-(1:3)], 1,
                                function(x) {Chat.Ind(x, sum(x))} )
   
   # extrapolated coverage with sample size * 2
   cov_extra <- apply(dat_wide[,-(1:3)], 1,
                      function(x) {Chat.Ind(x, 2*sum(x))})
   
   # base coverage 
   dat_biodiv$cov_base <- min(max(dat_biodiv$coverage, na.rm = T),
                              min(cov_extra, na.rm = T))
    
   # get mobr indices
   mob <- calc_biodiv(abund_mat = as.data.frame(dat_wide[,-(1:3)]),
                      groups = rep("frag", nrow(dat_wide)),
                      index = c("S", "S_n", "S_asymp", "S_PIE"), 
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
   dat_biodiv$S_cov <- rep(NA, nrow(dat_biodiv)) 
   
   dat_biodiv$S_chao <- mob$value[mob$index == "S_asymp"] 
   
   S_cov_std <- lapply(data.frame(t(dat_wide[,-(1:3)])),
                       function(x) try(estimateD(x, datatype = "abundance",
                                                 base = "coverage",
                                                 level = dat_biodiv$cov_base[1],
                                                 conf = NULL)))
   succeeded <- !sapply(S_cov_std, is.error)
   
   if (sum(succeeded) > 0){
      dat_biodiv$S_cov[succeeded] <- sapply(S_cov_std[succeeded],"[[","q = 0")
   }
   
   # standardized abundance and richness
   if (data_set$sample_design[1] == "standardized_fragment"){ 
      dat_biodiv$N_std <- dat_biodiv$N
      dat_biodiv$S_std_2 <- dat_biodiv$S_std_1 <- dat_biodiv$S_obs
   } else {
   # pooled and standardized plots   
      
      dat_biodiv$S_std_1 <- NA
      
      rel_sample_eff <- dat_biodiv$sample_eff/min(dat_biodiv$sample_eff)
      dat_biodiv$N_std <- dat_biodiv$N/rel_sample_eff
   
      for (i in 1:nrow(dat_biodiv)){
         if (dat_biodiv$N_std[i] >= n_thres)
            dat_biodiv$S_std_1[i] <- rarefaction(dat_wide[i,-(1:3)],
                                                 method = "indiv",
                                                 effort = round(dat_biodiv$N_std[i]))
      }
      
      if (data_set$sample_design[1] == "pooled"){
         dat_biodiv$S_std_2 <- dat_biodiv$S_std_1
      }
      
      if (data_set$sample_design[1] == "standardized_subsamples"){
         # get abundance by subsample
         dat_abund_sub <- data_set %>%
            group_by(frag_id, sample_id, species) %>%
            summarise(abundance = sum(abundance)) %>%
            ungroup()
         
         dat_sub <- dat_abund_sub %>%
            group_by(frag_id, sample_id) %>%
            summarise(S = n()) %>% ungroup
         
         dat_frag_mean <- dat_sub %>%
            group_by(frag_id) %>%
            summarise(S_std_2 = mean(S))
         
         dat_biodiv <- left_join(dat_biodiv, dat_frag_mean)
      }
   }
   
   # set indices to NA when there are no individuals
   empty_plots <- dat_biodiv$N == 0
   dat_biodiv$S_PIE[empty_plots] <- NA

   # set indices to NA when they are Inf or NaN
   dat_biodiv[,6:ncol(dat_biodiv)] <- lapply(dat_biodiv[,6:ncol(dat_biodiv)], Inf_to_NA)
   
   # set coverage standardized richness equal to observed
   # when observed coverage and base coverage equal 1
   coverage_eq_1 <- dat_biodiv$coverage == 1 & dat_biodiv$cov_base == 1
   dat_biodiv$S_cov[coverage_eq_1] <- dat_biodiv$S_obs[coverage_eq_1]
   
   # add fragment information
   dat_out <- left_join(dat_frag, dat_biodiv)

   # beta-diversity partitioning ---------------------------------------------
   
   # standardized sampling ----
   if (data_set$sample_design[1] == "standardized_fragment"){
      
      beta_div_comp <- get_beta_part(dat_wide) 
      
      beta_part_frag <- beta_div_comp$fragments 
      beta_part_study <- beta_div_comp$study 
   }
   
   # standardized subplots ----
   if (data_set$sample_design[1] == "standardized_subsamples"){
      
      # sum abundances in subsamples
      dat_abund_sub <- data_set %>%
         group_by(frag_id, frag_size_num, sample_id, species) %>%
         summarise(abundance = sum(abundance)) %>%
         ungroup()
      
      dat_wide_sub <- dat_abund_sub %>% 
         spread(key = species, value = abundance, fill = 0) 
      
      n_sub <- min(dat_wide_sub %>% 
                   group_by(frag_id) %>%
                   count() %>% 
                   ungroup() %>%
                   select(n)) 
      
      beta_part_samples <- map(1:n_resamples, betapart_subplots, dat1 = dat_wide_sub, n = n_sub)
      
      beta_part_frag <- beta_part_samples %>% map_dfr("fragments")
      beta_part_study <- beta_part_samples %>% map_dfr("study") 
   }
   
   # pooled samples ----
   if (data_set$sample_design[1] == "pooled"){
   
      beta_part_samples <- map(1:n_resamples, betapart_resample,
                               dat1 = dat_abund,
                               N_std = round(dat_biodiv$N_std))
      
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
   dat_frag2 <- select(dat_frag, frag_id, frag_size_char, frag_size_num)
   
   # dat_frag$site_label <- paste("Site", 1:nrow(dat_frag), sep = "")
   beta_part_frag <- beta_part_frag %>%
      left_join(dat_frag2, by = c("frag_x" = "frag_id")) %>%
      left_join(dat_frag2, by = c("frag_y" = "frag_id")) %>%
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
                    "betapart_study" = beta_part_study)
   
   return(out_list)
}    

# Execution of script -----------------------------------------------------

# read long format data file
infile <- path2Dropbox %+% "files_datapaper/Long_format_database/fragSAD_predicts_ewers.csv"
dat_long <- read.csv(infile, stringsAsFactors = F)
dim(dat_long)
str(dat_long)

head(dat_long)

# dat_long %>% select(dataset_label, sample_design) %>% distinct()
# 
# data_set <- dat_long %>% filter(dataset_label == "Ewers_2007")
# test <- get_biodiv(data_set)

parset <- expand.grid(fac_cont = c(2,10,100),
                      method_abund = c("as_is","round","ceiling","multiply"),
                      stringsAsFactors = F)

for (i in 1:nrow(parset)){
   out1 <- dat_long %>%
      split(.$dataset_label) %>%
      map(get_biodiv,
          fac_cont = parset$fac_cont[i],
          method_abund = parset$method_abund[i])
   
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
   
# out1 %>% 
#    select(dataset_label, sample_design) %>% 
#    distinct() %>%
#    ungroup() %>%
#    count(sample_design)
   
 
# base R version
# out1 <- by(dat_long, INDICES = list(dat_long$dataset_id)
#             FUN = get_biodiv)

# class(out1) <- "list"
# out1 <- bind_rows(out1)

# # purrr version
# out1 <- dat_long %>%
#    split(.$dataset_label) %>%
#    map_dfr(get_biodiv)
# 
# out1 %>% 
#    select(dataset_label, sample_design) %>% 
#    distinct() %>%
#    ungroup() %>%
#    count(sample_design)
# 
# # prepare output date
# path2outfile <- path2Dropbox %+% "files_datapaper/Analysis/biodiv_fragment_level.csv"
# write_csv(out1, path2outfile)

###
# check what happens with multiplication of abundances

# data_set1 <- dat_long %>% filter(dataset_label == "Cosson_1999")
# data_set3 <- data_set2 <- data_set1
# 
# data_set2$abundance <- 10*data_set2$abundance
# data_set3$abundance <- 100*data_set1$abundance
# 
# biodiv1 <- get_biodiv(data_set1)
# biodiv2 <- get_biodiv(data_set2)
# biodiv3 <- get_biodiv(data_set3)

# test <- out1 %>% filter(coverage == 1 & cov_base == 1)

