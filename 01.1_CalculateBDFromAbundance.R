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
# Function to calculate biodiversity indices from every patch

CalcBDfromAbundance <- function(filename){
   
   print(filename)
   
   filename2 <- strsplit(filename, split = "[.]")[[1]][1]
   
   dat_head <- read.xlsx(filename, sheetIndex = 1, rowIndex = 2:4,
                         header = F)
   dat_head_t <- as.data.frame(t(dat_head[,-1]))
   names(dat_head_t) <- dat_head[,1]
   
   if(filename %in% c("Dauber_2006_rounded.xls","Didham_1998.xls")){
      dat_head_t$entity.id <- sapply(as.character(dat_head_t$entity.id), function(x) strsplit(x,".",fixed=T)[[1]][1])
   }
   if(filename == "Schnitzler_2008_Tab2.xlsx"){
      dat_head_t$entity.id <- rep(paste0("fragment",1:10), times=c(rep(3,8),2,3))
   }
   
   
   
   ### extract fragment sizes and size ranks
   dat_frag_size <- read.xlsx(filename, sheetIndex = 1, rowIndex = 5,
                              stringsAsFactors = F, header = F)
   dat_frag_size <- dat_frag_size[, -1]
   
   dat_head_t$entity.size <- as.numeric(dat_frag_size[1,])
   
   dat_ranks <- read.xlsx(filename, sheetIndex = 1, rowIndex = 6, header = F,
                          stringsAsFactors = F)
   
   dat_ranks <- as.numeric(dat_ranks[1,-1])
   dat_ranks[is.na(dat_ranks)] <- 0
   
   dat_head_t$entity.size.rank <- dat_ranks
   
   ### extract sampling effort
   dat_sample_eff <- read.xlsx(filename, sheetIndex = 1, rowIndex = 1,
                               header = F, stringsAsFactors = F)
   dat_sample_eff <- as.numeric(dat_sample_eff[, -1])
   dat_sample_eff[is.na(dat_sample_eff)] <- 1
   
   if (max(dat_sample_eff) > 1) print(paste("Unequal sampling effort in", filename))
   
   ### extract abundance data
   dat_abund <- read.xlsx(filename, sheetIndex = 1, startRow = 7, header = F)
   dat_abund <- dat_abund[, -1]
   
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
   
   div_indi <- join(div_indi, dat_head_t[,c("entity.id","entity.size")], match="first")
   
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
   rel_sample_eff <- div_indi$sample_effort/max(div_indi$sample_effort)
   div_indi$N <- colSums(dat_abund_pool2)
#   div_indi$N_std <-  div_indi$N/div_indi$sample_effort ## standardize abundance by absolute sampling effort
   div_indi$N_std <-  div_indi$N/rel_sample_eff ## standardize abundance by relative sampling effort
   div_indi$S_obs <- colSums(dat_abund_pool2 > 0)
   S_plot <- colSums(dat_abund > 0)    # observed species richness per plot before plots within the same fragments have been pooled
   div_indi$S_plot <- ifelse(div_indi$sample_design == "pooled", NA, aggregate(S_plot, 
                            by = list(dat_head_t$entity.id,
                                      dat_head_t$entity.size.rank),
                            FUN = mean)$x) # mean of plot level species richness, comparison only meaningful if sampling effort among plots is equal, i.e. sampling_design either standardized or multiple standardized subsamples within fragments
   div_indi$ENS_pie <- diversity(t(dat_abund_pool2), index = "invsimpson")
   
   # sample coverage
   div_indi$coverage <- apply(dat_abund_pool2, 2,
                              function(x) {Chat.Ind(x, sum(x))} )
   
   # extrapolated coverage with sample size * 2
   div_indi$cov_extra <- apply(dat_abund_pool2, 2,
                      function(x) {Chat.Ind(x, 2*sum(x))})
   
   ### get base coverage following Chao et al. 2014. Ecol Monographs, box 1, p 60
   ### p.62:  However, for q = 0, extrapolation is reliable up to no more than double the reference sample size. Beyond that, the predictor for q = 0 may be subject to some bias because our asymptotic estimator for species richness (Chaol for abundance data and Chao2 for incidence data) is a lower bound only.
   div_indi$base_cov <- min(max(div_indi$coverage, na.rm = T),
                            min(div_indi$cov_extra, na.rm = T))

   # calculate standardized coverage
   div_indi$D0_hat <- rep(NA, nrow(div_indi)) 

   D_cov_std <- lapply(dat_abund_pool2,
                       function(x) try(estimateD(x, datatype = "abundance",
                                       base = "coverage",
                                       level = div_indi$base_cov[1],
                                       conf = NULL)))
   succeeded <- !sapply(D_cov_std, is.error)
   
   if (sum(succeeded) > 0){
      div_indi$D0_hat[succeeded] <- sapply(D_cov_std[succeeded],"[[","q = 0")
   }
   
   cov_eq_1 <- div_indi$coverage >= 1.0 
   cov_eq_1[is.na(cov_eq_1)] <- FALSE
   div_indi$D0_hat[cov_eq_1] <- div_indi$S_obs[cov_eq_1]

   # set indices to NA when there are no individuals
   empty_plots <- div_indi$N == 0 
   div_indi$ENS_pie[empty_plots] <- NA

   #############################################################################
   ### Beta-diversity partitioning
   ### Estimate species turnover that is due to replacement (in contrast to nestedness)
   
   div_indi$repl_part_BS_qT <- NA ##"BS" – Baselga family, Sørensen-based indices, computes quantitative forms of replacement, nestedness and D

   # Check sampling design
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
   
   if (div_indi$sample_design[1] != "pooled"){
      
      # Select fragments with 
      # Lower and upper quartile
      q_ranks <- quantile(dat_ranks, prob = c(0.25, 0.75))
      sub_ranks <- sort(dat_ranks[dat_ranks <= q_ranks[1] | dat_ranks >= q_ranks[2]])
      
      dat_abund_pool3 <- dat_abund_pool[dat_abund_pool[,2] %in% sub_ranks, ]
      
      size_class <- ifelse(dat_abund_pool3$Group.2 <= q_ranks[1], "Small", "Large")

      # divide by sampling effort in case of non-standardized design
      ncol1 <- ncol(dat_abund_pool3)
      if (div_indi$sample_design[1] == "subsamples_in_frag")
         dat_abund_pool3[,4:ncol1] <- dat_abund_pool3[,4:ncol1]/dat_abund_pool3$dat_sample_eff
            
      #Pool fragments with the same size rank
      dat_abund_pool3a <- aggregate(dat_abund_pool3[, 4:ncol1],
                                    by = list(size_class),
                                    FUN = mean)
      
      div_indi$repl_part_BS_qF <- beta.div.comp(dat_abund_pool3a[,-1], coef = "BS", quant = F)$part[4] # beta-div partitioning based on presence-absence, i.e. all species are equal regardless of their abundance
      div_indi$repl_part_BS_qT <- beta.div.comp(dat_abund_pool3a[,-1], coef = "BS", quant = T)$part[4] # beta-div partitioning based on abundance
   }
   
   # set indices to NA when they are Inf or NaN
   div_indi[,9:ncol(div_indi)] <- lapply(div_indi[,9:ncol(div_indi)], Inf_to_NA)
   
   return(div_indi)
}    

################################################################################
# execution of script

setwd(path2Dropbox %+% "good_datasets/")

div_list <- list()

filenames <- list.files(pattern="*.xls*", full.names = F)

filenames2 <- sapply(strsplit(filenames, split = "[.]"), "[[", 1)

for (i in 1:length(filenames)){
   temp <- try(CalcBDfromAbundance(filenames[i]))
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

