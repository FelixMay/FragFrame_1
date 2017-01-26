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
   
   # reading and cleaning the data
   dat_head <- read.xlsx(filename, sheetIndex = 1, startRow = 1, endRow = 3,
                         header = F)
   dat_head_t <- as.data.frame(t(dat_head[,-1]))
   names(dat_head_t) <- dat_head[,1]
   
   #dat_head <- dat_head[, names(dat_head) != "NA."]
   
   dat_frag_size <- read.xlsx(filename, sheetIndex = 1, startRow = 3, endRow = 4,
                              stringsAsFactors = F)
   dat_frag_size <- dat_frag_size[, -1]
   
   dat_head_t$entity.size <- as.numeric(dat_frag_size[1,])
   
   dat_ranks <- read.xlsx(filename, sheetIndex = 1, startRow = 4, endRow = 5)
   dat_ranks <- dat_ranks[, -1]
   #dat_ranks <- dat_ranks[, names(dat_ranks) != "NA."]
   
   dat_head_t$entity.size.rank <- as.numeric(dat_ranks[1,])
   
   # # just a quick check
   # plot(entity.size.rank ~ entity.size,
   #      data = dat_head_t[order(dat_head_t$entity.size), ], type = "b")

   dat_abund <- read.xlsx(filename, sheetIndex = 1, startRow = 6, header = F)
   dat_abund <- dat_abund[, -1]
   
   na_col <- apply(dat_abund, 1, function(x) sum(is.na(x)))
   na_row <- apply(dat_abund, 2, function(x) sum(is.na(x)))
   
   dat_abund <- dat_abund[na_col < dim(dat_abund)[2], na_row < dim(dat_abund)[1]]
   dat_abund[is.na(dat_abund)] <- 0
   
   # round abundances to integer numbers
   dat_abund <- round(dat_abund, digits = 0)
   dat_abund_t <- t(dat_abund)
   
   # pool data from the same fragments
   dat_abund_pool <- aggregate(dat_abund_t,
                               by = list(dat_head_t$entity.id,
                                         dat_head_t$entity.size.rank),
                               FUN = sum)
   dat_abund_pool2 <- as.data.frame(t(dat_abund_pool[,-c(1:2)]))
   names(dat_abund_pool2) <- dat_abund_pool[ ,1]
   
   # prepare output data
   div_indi <- data.frame(filename   = filename2, 
                          entity.id = dat_abund_pool[ ,1],
                          entity.size.rank =  dat_abund_pool[ ,2])
   rm(dat_abund_pool)
   
   # simple diversity indices
   div_indi$N <- colSums(dat_abund_pool2)
   div_indi$S <- colSums(dat_abund_pool2 > 0)
   
   div_indi$Shannon <- diversity(t(dat_abund_pool2), index = "shannon")
   div_indi$PIE <- diversity(t(dat_abund_pool2), index = "simpson")
   
   div_indi$ENS_shannon <- exp(div_indi$Shannon)
   div_indi$ENS_pie <- diversity(t(dat_abund_pool2), index = "invsimpson")
   
   div_indi$Pielou_even <- div_indi$Shannon / log(div_indi$S)
   
   # sample coverage
   div_indi$coverage <- apply(dat_abund_pool2, 2,
                              function(x) {Chat.Ind(x, sum(x))} )
   
   # extrapolated coverage with sample size * 2
   cov_extra <- apply(dat_abund_pool2, 2,
                      function(x) {Chat.Ind(x, 2*sum(x))})
   
   # get base coverage following Chao et al. 2014. Ecol Monographs
   div_indi$base_cov <- max(max(div_indi$coverage, na.rm = T),
                            min(cov_extra, na.rm = T))
   
   # get more conservative base coverage - at maximum extrapolate to 2*n
   # div_indi$base_cov <- min(max(div_indi$coverage, na.rm = T),
   #                          min(cov_extra, na.rm = T))
   
   # calculate standardized coverage
   div_indi$D0_hat <- rep(NA, nrow(div_indi)) 
   div_indi$D1_hat <- rep(NA, nrow(div_indi)) 
   div_indi$D2_hat <- rep(NA, nrow(div_indi)) 
   
   D_cov_std <- lapply(dat_abund_pool2,
                       function(x) try(estimateD(x, datatype = "abundance",
                                       base = "coverage",
                                       level = div_indi$base_cov[1],
                                       conf = NULL)))
   succeeded <- !sapply(D_cov_std, is.error)
   
   if (sum(succeeded) > 0){
      div_indi$D0_hat[succeeded] <- sapply(D_cov_std[succeeded],"[[","q = 0")
      div_indi$D1_hat[succeeded] <- sapply(D_cov_std[succeeded],"[[","q = 1")
      div_indi$D2_hat[succeeded] <- sapply(D_cov_std[succeeded],"[[","q = 2")
   }
   
   cov_eq_1 <- div_indi$coverage >= 1.0 
   cov_eq_1[is.na(cov_eq_1)] <- FALSE
   div_indi$D0_hat[cov_eq_1] <- div_indi$S[cov_eq_1]
   div_indi$D1_hat[cov_eq_1] <- div_indi$ENS_shannon[cov_eq_1]
   div_indi$D2_hat[cov_eq_1] <- div_indi$ENS_pie[cov_eq_1]
   
   # set indices to NA when there are no individuals
   empty_plots <- div_indi$N == 0 
   div_indi$Shannon[empty_plots] <- NA
   div_indi$PIE[empty_plots] <- NA
   div_indi$ENS_shannon[empty_plots] <- NA
   div_indi$ENS_pie[empty_plots] <- NA
   div_indi$Pielou_even[empty_plots] <- NA
   
   # extrapolation to asymptotic species richness
   D_asymp_list <- lapply(dat_abund_pool2,
                          function(x) try(SpadeR::Diversity(x, datatype = "abundance",
                                                            q = c(0,1,2))))
   succeeded <- !sapply(D_asymp_list, is.error)
   D_asymp_mat <- matrix(NA, nrow = ncol(dat_abund_pool2), ncol = 8)
   colnames(D_asymp_mat) <- c("Chao1", "Chao1-bc", "iChao1", "ACE", "ACE-1",
                              "D0_asymp","D1_asymp","D2_asymp")   

   if (sum(succeeded) > 0){
      S_asymp <- sapply(D_asymp_list[succeeded], function(div1){div1$SPECIES.RICHNESS[,"Estimate"]})
      D_asymp_mat[succeeded, 1:5] <- t(S_asymp)
      
      Hill_asymp <- sapply(D_asymp_list[succeeded], function(div1){div1$HILL.NUMBERS[,"ChaoJost"]})
      D_asymp_mat[succeeded, 6:8] <- t(Hill_asymp)
   }

   div_indi <- cbind(div_indi, D_asymp_mat)
   
   # # inter- and extrapolation plot
   # inext1 <- iNEXT(dat_abund_pool2[ ,div_indi$N > 0], q = 0, datatype="abundance")
   # plot1 <- ggiNEXT(inext1, type = 1)
   # #plot2 <- ggiNEXT(inext1, type = 3)
   # 
   # # save plot and summary statistics
   # fig_name <- paste(path2temp, filename2, ".pdf", sep="")
   # pdf(fig_name, width = 7, height = 7)
   # #grid.arrange(plot1, plot2, ncol = 2)
   # print(plot1)
   # dev.off()
   
   # set indices to NA when they are Inf or NaN
   div_indi[,9:23] <- lapply(div_indi[,9:23], Inf_to_NA)
   
   return(div_indi)
}    

################################################################################
# execution of script

setwd(path2Dropbox %+% "good_datasets/")

div_list <- list()

filenames <- list.files(pattern="*.xls*", full.names = F)

# remove data set with non-integer numbers
# filenames <- filenames[filenames != "De_Lima_1999.xlsx" &
#                        filenames != "Leal_2012.xls" & 
#                        filenames != "Gavish_2012.xlsx" &
#                        filenames != "Gavish_2012_B.xlsx"]
filenames2 <- sapply(strsplit(filenames, split = "[.]"), "[[", 1)


# try if all files really work
for (i in 1:length(filenames)){
   temp <- CalcBDfromAbundance(filenames[i])
   div_list[[filenames2[i]]] <- temp
}

# for (i in 1:length(filenames)){
#    temp <- try(CalcBDfromAbundance(filenames[i]))
#    if (!inherits(temp, "try-error")){
#       div_list[[filenames2[i]]] <- temp
#    }
# }   

div_df <- bind_rows(div_list)

### get rid of irrelevant data, e.g. matrix, clearcut
div_df_nomatrix <- filter(div_df, entity.size.rank > 0)

write.table(div_df_nomatrix, file = paste(path2temp, "DiversityData.csv", sep = ""),
            sep = ";", row.names = F)

# check incomplete cases
summary(div_df_nomatrix)
div_df_nomatrix[!complete.cases(div_df_nomatrix), ]

div_df_nomatrix[is.na(div_df_nomatrix$Pielou_even), ]
div_df_nomatrix[is.na(div_df_nomatrix$D1_hat), ]
