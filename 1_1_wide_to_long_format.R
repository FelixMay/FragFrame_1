
############################################
# Read files, transpose and switch to long format

read_data_files <- function(filename){
   
   print(filename)
   
   # filename2 <- strsplit(filename, split = "[.]")[[1]][1]

   path2file <- path2Dropbox %+% "files_datapaper/" %+% filename
   
   # Sampling effort
   sample_eff <- read.table(path2file, sep = ",", row.names = 1, nrows = 1)
   sample_eff <- sample_eff[-length(sample_eff)]
   sample_eff <- as.numeric(sample_eff)
   
   if (max(sample_eff) > 1) print(paste("Unequal sampling effort in", filename))
   
   # entity_id_orig
   entity_id_orig <- read.table(path2file, sep = ",", row.names = 1, nrows = 1,
                                skip = 1, stringsAsFactors = F)
   entity_id_orig <- as.character(entity_id_orig[1,1:length(sample_eff)])
   
   # entity_id_plot
   entity_id_plot <- read.table(path2file, sep = ",", row.names = 1, nrows = 1,
                                skip = 2, stringsAsFactors = F)
   entity_id_plot <- as.numeric(entity_id_plot[1,1:length(sample_eff)])
   
   # entity_size
   entity_size <- read.table(path2file, sep = ",", row.names = 1, nrows = 1,
                             skip = 3, stringsAsFactors = F)
   entity_size <- as.numeric(entity_size[1,1:length(sample_eff)])
   
   # abundance data
   dat_abund <- read.table(path2file, sep = ",", row.names = 1,
                           skip = 4, stringsAsFactors = F)
   dat_abund <- dat_abund[,1:length(sample_eff)]
   dat_abund_t <- as.data.frame(t(dat_abund))
   
   # combine fragment-data
   dat1 <- data.frame(entity_id_orig, entity_id_plot, entity_size, sample_eff)
   
   # determine sampling design
   sample_eff_per_frag <- tapply(dat1$sample_eff, dat1$entity_id_orig, sum)
   range_sample_eff <- range(sample_eff_per_frag)   
   range_sample_units <- range(entity_id_plot)
   
   if (range_sample_eff[2] - range_sample_eff[1] == 0){
      dat1$sample_design <- "standardized"
   } else {
      if (range_sample_units[2] - range_sample_units[1] == 0) {
         dat1$sample_design <- "pooled"
      } else {
         dat1$sample_design <- "subsamples_in_frag"
      }
   }
   
   
   
   
   
   ### extract abundance data
   dat_abund <- read.xlsx(path2file, sheetIndex = 1, startRow = 7, header = F)
   dat_abund <- as.data.frame(dat_abund[, -1])
   
   na_col <- apply(dat_abund, 1, function(x) sum(is.na(x)))
   na_row <- apply(dat_abund, 2, function(x) sum(is.na(x)))
   
   dat_abund <- dat_abund[na_col < dim(dat_abund)[2], na_row < dim(dat_abund)[1]]
   dat_abund[is.na(dat_abund)] <- 0
   
   # round abundances to nearest upper integer
   dat_abund <- ceiling(dat_abund)
   dat_abund_t <- t(dat_abund)
   
   dat_abund_t <- cbind(dat_sample_eff, dat_abund_t)
   
   ### pool data from the same fragments
   dat_abund_pool <- aggregate(dat_abund_t,
                               by = list(dat_head_t$entity.id,
                                         dat_head_t$entity.size.rank),
                               FUN = sum)
   dat_abund_pool2 <- as.data.frame(t(dat_abund_pool[,-c(1:3)]))
   names(dat_abund_pool2) <- dat_abund_pool[ ,1]
   
   ### get number of samples per fragments
   dat_n_sampling_units <- aggregate(dat_head_t$entity.id.plot,
                                     by = list(dat_head_t$entity.id,
                                               dat_head_t$entity.size.rank),
                                     FUN = length)
   
   ### prepare output data
   div_indi <- data.frame(filename   = filename2, 
                          entity.id = dat_abund_pool[ ,1],
                          entity.size.rank =  dat_abund_pool[ ,2])
   
   div_indi <- div_indi %>% 
      left_join(dat_head_t[,c("entity.id","entity.size","entity.type")]) %>%
      distinct()
   
   ### Check sampling design
   div_indi$sampling_units <- dat_n_sampling_units[,3]
   div_indi$sample_effort <- dat_abund_pool[,"dat_sample_eff"] 
   
   range_sample_eff <- range(div_indi$sample_effort)
   range_sample_units <- range(div_indi$sampling_units)
   if (range_sample_eff[2] - range_sample_eff[1] == 0){
      div_indi$sample_design <- "standardized"
   } else {
      if (range_sample_units[2] - range_sample_units[1] == 0) {
         div_indi$sample_design <- "pooled"
      } else {
         div_indi$sample_design <- "subsamples_in_frag"
      }
   }
   
   ### simple diversity indices
   
   # standardize sampling effort
   # smallest sampling unit:  sampling effort = 1
   rel_sample_eff <- div_indi$sample_effort/min(div_indi$sample_effort)
   
   # No. of individuals
   div_indi$N <- colSums(dat_abund_pool2)
   
   # n_min <- min(div_indi$N)
   
   # No. of individuals standardized by sampling effort
   
   #div_indi$N_std <-  div_indi$N/div_indi$sample_effort ## standardize abundance by absolute sampling effort
   div_indi$N_std <-  div_indi$N/rel_sample_eff ## standardize abundance by relative sampling effort
   
   # Determine base sample size for the sample-size-based rarefaction and extrapolation
   # Chao et al. 2014. Box 1
   n_max <- max(div_indi$N)
   n_min_r <- min(2 * div_indi$N)
   
   n_base <- min(n_max, n_min_r)   
   div_indi$n_base <- n_base
   
   # Chao et al. 2014 suggests using the maximum here,
   # we decided to go for a more conservative measure
   
   # Determine base coverage for the coverage-based rarefaction and extrapolation
   # Chao et al. 2014. Box 1
   
   # sample coverage
   div_indi$coverage <- apply(dat_abund_pool2, 2,
                              function(x) {Chat.Ind(x, sum(x))} )
   
   # extrapolated coverage with sample size * 2
   cov_extra <- apply(dat_abund_pool2, 2,
                      function(x) {Chat.Ind(x, 2*sum(x))})
 
   # base coverage 
   cov_base <- min(max(div_indi$coverage, na.rm = T),
                   min(cov_extra, na.rm = T))
   div_indi$cov_base <- cov_base
   # get mobr indices
   div_dat <- calc_biodiv(abund_mat = t(dat_abund_pool2),
                          groups = rep("frag", ncol(dat_abund_pool2)),
                          index = c("N", "S", "S_n", "S_asymp", "f_0", "S_PIE"), 
                          effort = n_base,
                          extrapolate = T,
                          return_NA = F)
   
   # observed richness
   div_indi$S_obs <- div_dat$value[div_dat$index == "S"]
   
   # richness standardized by mean sampling effort in smallest sampling unit
   div_indi$S_std <- NA
   for (i in 1:nrow(div_indi)){
      if (div_indi$N_std[i] >= n_thres)   
         div_indi$S_std[i] <- rarefaction(dat_abund_pool2[,i],
                                          method = "indiv",
                                          effort = round(div_indi$N_std[i]))        
   }

   # rarefied richness
   if (n_base >= n_thres){
      div_indi$S_n <- div_dat$value[div_dat$index == "S_n"]
   } else {
      div_indi$S_n <- NA
   }
   
   # asymptotic richness
   div_indi$S_asymp <- div_dat$value[div_dat$index == "S_asymp"]
   
   # ENS PIE
   div_indi$S_PIE <- div_dat$value[div_dat$index == "S_PIE"]
   
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
   
   # Calculate coverage standardized richness
   div_indi$S_cov <- rep(NA, nrow(div_indi)) 
   
   S_cov_std <- lapply(dat_abund_pool2,
                       function(x) try(estimateD(x, datatype = "abundance",
                                       base = "coverage",
                                       level = cov_base,
                                       conf = NULL)))
   succeeded <- !sapply(S_cov_std, is.error)

   if (sum(succeeded) > 0){
      div_indi$S_cov[succeeded] <- sapply(S_cov_std[succeeded],"[[","q = 0")
   }
   
   # I am not sure these three lines are neede?
   # cov_eq_1 <- div_indi$coverage >= 1.0
   # cov_eq_1[is.na(cov_eq_1)] <- FALSE
   # div_indi$D0_hat[cov_eq_1] <- div_indi$S_obs[cov_eq_1]

   # set indices to NA when there are no individuals
   empty_plots <- div_indi$N == 0
   div_indi$S_asymp[empty_plots] <- NA
   div_indi$S_PIE[empty_plots] <- NA

   #############################################################################
   ### Beta-diversity partitioning
   ### Estimate species turnover that is due to replacement (in contrast to nestedness)
   
   J_qF <- beta.div.comp(t(dat_abund_pool2), coef = "J", quant = F)
   S_qF <- beta.div.comp(t(dat_abund_pool2), coef = "S", quant = F)
   BJ_qF <- beta.div.comp(t(dat_abund_pool2), coef = "BJ", quant = F)
   BS_qF <- beta.div.comp(t(dat_abund_pool2), coef = "BS", quant = F)
   
   J_qT <- beta.div.comp(t(dat_abund_pool2), coef = "J", quant = T)
   S_qT <- beta.div.comp(t(dat_abund_pool2), coef = "S", quant = T)
   BJ_qT <- beta.div.comp(t(dat_abund_pool2), coef = "BJ", quant = T)
   BS_qT <- beta.div.comp(t(dat_abund_pool2), coef = "BS", quant = T)
   
   beta_div_tab <- dist_to_dataframe(J_qF$repl)
   names(beta_div_tab)[3] <- "J_qF_repl"
   
   beta_div_tab$J_qF_rich <- J_qF$rich
   
   beta_div_tab$S_qF_repl <- S_qF$repl
   beta_div_tab$S_qF_rich <- S_qF$rich
   
   beta_div_tab$BJ_qF_repl <- BJ_qF$repl
   beta_div_tab$BJ_qF_rich <- BJ_qF$rich
   
   beta_div_tab$BS_qF_repl <- BS_qF$repl
   beta_div_tab$BS_qF_rich <- BS_qF$rich
   
   beta_div_tab$J_qT_repl <- J_qF$repl
   beta_div_tab$J_qT_rich <- J_qF$rich
   
   beta_div_tab$S_qT_repl <- S_qT$repl
   beta_div_tab$S_qT_rich <- S_qT$rich
   
   beta_div_tab$BJ_qT_repl <- BJ_qT$repl
   beta_div_tab$BJ_qT_rich <- BJ_qT$rich
   
   beta_div_tab$BS_qT_repl <- BS_qT$repl
   beta_div_tab$BS_qT_rich <- BS_qT$rich
   
   
   
   # # old version
   # div_indi$repl_part_BS_qF <- NA # "BS" – Baselga family, Sørensen-based indices, computes presence/absence form
   # div_indi$repl_part_BS_qT <- NA # computes quantitative form
   # 
   # # Check sampling design
   # range_sample_eff <- range(div_indi$sample_effort)
   # range_sample_units <- range(div_indi$sampling_units)
   # if (range_sample_eff[2] - range_sample_eff[1] == 0){
   #    div_indi$sample_design <- "standardized"
   # } else {
   #    if (range_sample_units[2] - range_sample_units[1] == 0) {
   #       div_indi$sample_design <- "pooled"
   #    } else {
   #       div_indi$sample_design <- "subsamples_in_frag"
   #    }
   # }
   # 
   # if (div_indi$sample_design[1] != "pooled"){
   #    
   #    # Select fragments with 
   #    # Lower and upper quartile
   #    q_ranks <- quantile(dat_ranks, prob = c(0.25, 0.75))
   #    sub_ranks <- sort(dat_ranks[dat_ranks <= q_ranks[1] | dat_ranks >= q_ranks[2]])
   #    
   #    dat_abund_pool3 <- dat_abund_pool[dat_abund_pool[,2] %in% sub_ranks, ]
   #    
   #    size_class <- ifelse(dat_abund_pool3$Group.2 <= q_ranks[1], "Small", "Large")
   # 
   #    # divide by sampling effort in case of non-standardized design
   #    ncol1 <- ncol(dat_abund_pool3)
   #    if (div_indi$sample_design[1] == "subsamples_in_frag")
   #       dat_abund_pool3[,4:ncol1] <- dat_abund_pool3[,4:ncol1]/dat_abund_pool3$dat_sample_eff
   #          
   #    #Pool fragments with the same size rank
   #    dat_abund_pool3a <- aggregate(dat_abund_pool3[, 4:ncol1],
   #                                  by = list(size_class),
   #                                  FUN = mean)
   #    
   #    div_indi$repl_part_BS_qF <- beta.div.comp(dat_abund_pool3a[,-1], coef = "BS", quant = F)$part[4] # beta-div partitioning based on presence-absence, i.e. all species are equal regardless of their abundance
   #    div_indi$repl_part_BS_qT <- beta.div.comp(dat_abund_pool3a[,-1], coef = "BS", quant = T)$part[4] # beta-div partitioning based on abundance
   # }
   
   # set indices to NA when they are Inf or NaN
   div_indi[,9:ncol(div_indi)] <- lapply(div_indi[,9:ncol(div_indi)], Inf_to_NA)
   
   return(div_indi)
}    

################################################################################
# execution of script

div_list <- list()

filenames <- list.files(path =  path2Dropbox %+% "files_datapaper/",
                        pattern="*.csv", full.names = F)

# filenames2 <- sapply(strsplit(filenames, split = "[.]"), "[[", 1)

for (i in 1:length(filenames)){
   temp <- try(CalcBDfromAbundance(filenames[i], n_thres = 5))
   if (!inherits(temp, "try-error")){
      div_list[[filenames2[i]]] <- temp
   }
}

div_df <- bind_rows(div_list)

### get rid of irrelevant data, e.g. matrix, clearcut
div_df_nomatrix <- filter(div_df, entity.size.rank > 0)

write.table(div_df_nomatrix, file = paste(path2temp, "DiversityData.csv", sep = ""),
            sep = ",", row.names = F)

# # check incomplete cases
# summary(div_df_nomatrix)
# div_df_nomatrix[!complete.cases(div_df_nomatrix), ]

