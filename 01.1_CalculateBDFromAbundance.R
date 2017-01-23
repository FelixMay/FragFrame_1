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
   
   apply(dat_abund_pool2, 2, function(x) {Chat.Ind(x, sum(x))} )
   div_indi$Sample_coverage <- apply(dat_abund_pool2, 2,
                                     function(x) {Chat.Ind(x, sum(x))} )
   
   # set indices to NA when there are no individuals
   empty_plots <- div_indi$N == 0 
   div_indi$Shannon[empty_plots] <- NA
   div_indi$PIE[empty_plots] <- NA
   div_indi$ENS_shannon[empty_plots] <- NA
   div_indi$ENS_pie[empty_plots] <- NA
   div_indi$Pielou_even[empty_plots] <- NA

   # # coverage-based indices
   # temp <- estimateD(dat_abund_pool2, "abundance", base="coverage", level=0.99, conf=NULL) # species richness, shannon, simpson standardized by coverage (99 %)
   # div_indi$S_cov <- temp[,"q = 0"]
   # div_indi$Shannon_cov <- temp[,"q = 1"]
   # div_indi$PIE_cov <- temp[,"q = 2"]

   # # rarefaction curves (=non-spatial accumulation curves )
   # #SAC_list <- lapply(dat_abund, SAC.coleman)
   # 
   # # rarefied richness at minimum number of individuals
   # n_min <- 10
   # n_min_sample = max(n_min, min(div_indi$N))
   # plots_low_n = div_indi$N < n_min
   # 
   # 
   # nmin <- min(div_indi$N)
   # #div_indi$S_rare1 <- sapply(SAC_list, "[", nmin)
   # div_indi$S_rare <- rarefy(dat_abund_pool2, nmin, MARGIN = 2)
   
   # # extrapolation
   # chao_list <- lapply(dat_abund_pool2, function(x) try(SpadeR::ChaoSpecies(x, datatype = "abundance")))
   # 
   # succeeded <- !sapply(chao_list, is.error)
   # chao_mat <- matrix(NA, nrow = 9, ncol = ncol(dat_abund_pool2))
   # rownames(chao_mat) <- c("Homogeneous_Model","Homogeneous_MLE", "Chao1",
   #                         "Chao1-bc", "iChao1", "ACE", "ACE-1" ,
   #                         "1st_order_jackknife","2nd_order_jackknife")   
   # 
   # if (sum(succeeded) > 0){
   #    chao_spec <- sapply(chao_list[succeeded], function(chao1){chao1$Species_table[,"Estimate"]})
   #    
   #    chao_mat[, succeeded] <- chao_spec
   # }
   # 
   # div_indi <- cbind(div_indi, t(chao_mat))
   
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


for (i in 1:length(filenames)){
   temp <- try(CalcBDfromAbundance(filenames[i]))
   if (!inherits(temp, "try-error")){
      div_list[[filenames2[i]]] <- temp
   }
}   
#div_list[[filenames2[i]]] <- div_indi

div_df <- bind_rows(div_list)

### get rid of irrelevant data, e.g. matrix, clearcut
div_df_nomatrix <- filter(div_df, entity.size.rank > 0)

write.table(div_df_nomatrix, file = paste(path2temp, "DiversityData.csv", sep = ""),
            sep = ";", row.names = F)

#write.csv(div_df_nomatrix, file = paste(path2temp, "DiversityData.csv", sep = ""))
            