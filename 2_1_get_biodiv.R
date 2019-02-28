# source code from iNEXT package
#
# Abundance-based sample coverage -----------------------------------------
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


# Sample subplots from community matrix
betapart_subplots <- function(index, dat1, n){
   dat_sample <- dat1 %>%
      sample_n(n)
   class(dat_sample) <- "data.frame"
   rownames(dat_sample) <- dat_sample$frag_id
   
   pars <- expand.grid(coef = c("S","J","BS","BJ"),
                       quant = c(FALSE,TRUE), stringsAsFactors = F)
   
   beta_part_list <- list()
   for (i in 1:nrow(pars)){
      beta_part <- beta.div.comp(dat_sample[,-(1:2)],
                                 coef = pars$coef[i],
                                 quant = pars$quant[i])
      beta_div_repl <- dist_to_dataframe(beta_part$repl)
      beta_div_rich <- dist_to_dataframe(beta_part$rich)
      beta_div_repl$part <- "repl"
      beta_div_rich$part <- "rich"
      beta_div <- rbind(beta_div_repl, beta_div_rich)
      beta_div$coef <- pars$coef[i]
      beta_div$quant <- pars$quant[i]
      
      beta_part_list[[i]] <- beta_div
   }
}

# Function to calculate biodiversity indices for every fragment  ----------
get_biodiv <- function(data_set, n_thres = 5, fac_cont = 10, method_abund = c("as_is","round","ceiling","multiply")){
   
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
   
   dat_wide <- dat_frag %>%
      select(frag_id, sample_eff) %>%
      left_join(dat_wide)
   
   
   # get biodiversity indices ----
   dat_biodiv <- data.frame(frag_id = dat_wide$frag_id,
                            sample_eff = dat_wide$sample_eff,
                            N = rowSums(dat_wide[,-(1:2)]),
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
   dat_biodiv$coverage <- apply(dat_wide[,-(1:2)], 1,
                                function(x) {Chat.Ind(x, sum(x))} )
   
   # extrapolated coverage with sample size * 2
   cov_extra <- apply(dat_wide[,-(1:2)], 1,
                      function(x) {Chat.Ind(x, 2*sum(x))})
   
   # base coverage 
   dat_biodiv$cov_base <- min(max(dat_biodiv$coverage, na.rm = T),
                              min(cov_extra, na.rm = T))
    
   # get mobr indices
   mob <- calc_biodiv(abund_mat = as.data.frame(dat_wide[,-(1:2)]),
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
   dat_biodiv$S_cov <- rep(NA, nrow(dat_biodiv)) 
   
   S_cov_std <- lapply(data.frame(t(dat_wide[,-(1:2)])),
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
      rel_sample_eff <- dat_biodiv$sample_eff/min(dat_biodiv$sample_eff)
      dat_biodiv$N_std <- dat_biodiv$N/rel_sample_eff
   
      for (i in 1:nrow(dat_biodiv)){
         if (dat_biodiv$N_std[i] >= n_thres)
            dat_biodiv$S_std_1[i] <- rarefaction(dat_wide[i,-(1:2)],
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
   
   if (data_set$sample_design[1] == "standardized"){
      
      # add rownames for betadivpart
      class(dat_wide) <- "data.frame"
      rownames(dat_wide) <- dat_wide$frag_id
      
      J_qF <- beta.div.comp(dat_wide[,-(1:2)], coef = "J", quant = F)
      
      beta_div_tab <- dist_to_dataframe(J_qF$repl)
      names(beta_div_tab)[1:3] <- c("frag_x", "frag_y", "J_qF_repl")
      
      # beta_part_study <- data.frame(dataset_label = data_set$dataset_label[1],
      #                               repl_part = J_qF$part["Repl/BDtotal"],
      #                               richdiff_part = J_qF$part["RichDif/BDtotal"],
      #                               method = J_qF$Note,
      #                               stringsAsFactors = F, row.names = NULL)
   }
   
   if (data_set$sample_design[1] == "standardized_subsamples"){
      
      # sum abundances in subsamples
      dat_abund_sub <- data_set %>%
         group_by(frag_id, sample_id, species) %>%
         summarise(abundance = sum(abundance))
      
      dat_wide_sub <- dat_abund_sub %>% 
         spread(key = species, value = abundance, fill = 0) %>%
         group_by(frag_id)
      
      n_sub <- min(dat_wide_sub %>% 
                      count() %>% 
                      ungroup() %>%
                      select(n)) 
      
      beta_part_sample <- map_dfr(1:100, betapart_subplots, dat1 = dat_wide_sub, n = n_sub)
      
      beta_part_sample %>% group_by(row, col, part) %>%
         summarise(mean = mean(value))
      
     
      
      beta_div_tab <- dist_to_dataframe(J_qF$repl)
      names(beta_div_tab)[1:3] <- c("frag_x", "frag_y", "J_qF_repl")
   }
   
   dat_frag2 <- select(dat_frag, frag_id, frag_size_char, frag_size_num)
   
   # dat_frag$site_label <- paste("Site", 1:nrow(dat_frag), sep = "")
   beta_div_tab <- beta_div_tab %>%
      left_join(dat_frag2, by = c("frag_x" = "frag_id")) %>%
      left_join(dat_frag2, by = c("frag_y" = "frag_id")) %>%
      mutate(diff_area = frag_size_num.y - frag_size_num.x,
             log10_ratio_area = log10(frag_size_num.y/frag_size_num.x)
      )
   
   # adjust dataframe
   beta_div_tab$dataset_id <- data_set$dataset_id[1]
   beta_div_tab$dataset_label <- data_set$dataset_label[1]
   beta_div_tab$sample_design <- data_set$sample_design[1]
   
   beta_div_tab <- beta_div_tab %>% 
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
   
   out_list <- list("biodiv" = dat_out,
                    "betapart" = beta_div_tab)
   
   return(out_list)
}    

# execution of script -----------------------------------------------------

# read long format data file
infile <- path2Dropbox %+% "files_datapaper/Long_format_database/fragSAD_long.csv"
dat_long <- read.csv(infile, stringsAsFactors = F)
dim(dat_long)
str(dat_long)

head(dat_long)


data_set <- dat_long %>% filter(dataset_label == "Nogueira_2016")

get_biodiv(data_set)

parset <- expand.grid(fac_cont = c(2,10,100),
                      method_abund = c("as_is","round","ceiling","multiply"),
                      stringsAsFactors = F)

for (i in 1:nrow(parset)){
   out1 <- dat_long %>%
      split(.$dataset_label) %>%
      map(get_biodiv,
          fac_cont = parset$fac_cont[i],
          method_abund = parset$method_abund[i])
   
   out_biodiv <- out1 %>% map_dfr("biodiv")
   out_betapart <- out1 %>% map_dfr("betapart")
   
   # prepare output date
   outfile_name <- i %+% "_biodiv_frag_fcont_" %+% parset$fac_cont[i] %+%
      "_mabund_"%+% parset$method_abund[i] %+% ".csv"
   path2outfile <- path2Dropbox %+% "files_datapaper/Analysis/" %+% outfile_name
   write_csv(out_biodiv, path2outfile)
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

