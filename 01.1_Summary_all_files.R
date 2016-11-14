setwd(path2Dropbox %+% "good_datasets/")

div_list <- list()

filenames <- list.files(pattern="*.xls*", full.names = F)
filenames2 <- sapply(strsplit(filenames, split = "[.]"), "[[", 1)

# remove data set with non-integer numbers
filenames <- filenames[filenames != "Leal_2012_new_.xls"]

for (i in 1:length(filenames)){
   
   print(filenames[i])
   
   # reading and cleaning the data
   dat_head <- read.xlsx(filenames[i], sheetIndex = 1, startRow = 1, endRow = 3,
                         header = F)
   dat_head_t <- as.data.frame(t(dat_head[,-1]))
   names(dat_head_t) <- dat_head[,1]
   
   #dat_head <- dat_head[, names(dat_head) != "NA."]
   
   dat_ranks <- read.xlsx(filenames[i], sheetIndex = 1, startRow = 4, endRow = 5)
   dat_ranks <- dat_ranks[, -1]
   #dat_ranks <- dat_ranks[, names(dat_ranks) != "NA."]
   
   dat_head_t$entity.size.rank <- as.numeric(dat_ranks[1,])
   
   dat_abund <- read.xlsx(filenames[i], sheetIndex = 1, startRow = 6, header = F)
   dat_abund <- dat_abund[, -1]
   
   na_col <- apply(dat_abund, 1, function(x) sum(is.na(x)))
   na_row <- apply(dat_abund, 2, function(x) sum(is.na(x)))
   
   dat_abund <- dat_abund[na_col < dim(dat_abund)[2], na_row < dim(dat_abund)[1]]
   dat_abund[is.na(dat_abund)] <- 0
   
   dat_abund_t <- t(dat_abund)
   
   # pool data from the same fragments
   dat_abund_pool <- aggregate(dat_abund_t, by = list(dat_head_t$entity.id), FUN = sum)
   dat_abund_pool2 <- as.data.frame(t(dat_abund_pool[,-1]))
   names(dat_abund_pool2) <- dat_abund_pool[,1]
   rm(dat_abund_pool)
   
   # prepare output data
   dat_head_unique <- unique(dat_head_t[,c(1,4)])
   
   div_indi <- data.frame(filename   = filenames2[i], 
                          entity.id = dat_head_unique$entity.id,
                          entity.size.rank = dat_head_unique$entity.size.rank)
   
   # simple diversity indices
   div_indi$N <- colSums(dat_abund_pool2)
   div_indi$S <- colSums(dat_abund_pool2 > 0)
   
   div_indi$Shannon <- diversity(t(dat_abund_pool2), index = "shannon")
   div_indi$PIE <- diversity(t(dat_abund_pool2), index = "simpson")
   
   div_indi$ENS_shannon <- exp(div_indi$Shannon)
   div_indi$ENS_pie <- diversity(t(dat_abund_pool2), index = "invsimpson")
   
   # rarefaction curves (=non-spatial accumulation curves )
   #SAC_list <- lapply(dat_abund, SAC.coleman)
   
   # rarefied richness at minimum number of individuals
   nmin <- min(div_indi$N)
   #div_indi$S_rare1 <- sapply(SAC_list, "[", nmin)
   div_indi$S_rare <- rarefy(dat_abund_pool2, nmin, MARGIN = 2)
   
   # extrapolation
   chao_list <- lapply(dat_abund_pool2, function(x) try(SpadeR::ChaoSpecies(x, datatype = "abundance")))
   succeeded <- !sapply(chao_list, is.error)
   
   chao_mat <- matrix(NA, nrow = 9, ncol = ncol(dat_abund_pool2))
   chao_spec <- sapply(chao_list[succeeded], function(chao1){chao1$Species.Table[,"Estimate"]})
   
   chao_mat[, succeeded] <- chao_spec
   rownames(chao_mat) <- rownames(chao_spec)
   
   div_indi <- cbind(div_indi, t(chao_mat))
   
   div_list[[filenames2[i]]] <- div_indi
   
   # inter- and extrapolation plot
   inext1 <- iNEXT(dat_abund_pool2, q = 0, datatype="abundance")
   plot1 <- ggiNEXT(inext1, type = 1)
   #plot2 <- ggiNEXT(inext1, type = 3)
   
   # save plot and summary statistics
   fig_name <- paste(path2temp, filenames2[i], ".pdf", sep="")
   pdf(fig_name, width = 7, height = 7)
   #grid.arrange(plot1, plot2, ncol = 2)
   print(plot1)
   dev.off()
}

div_df <- bind_rows(div_list)
div_df_nomatrix <- filter(div_df, entity.size.rank > 0)

write.table(div_df, file = paste(path2temp, "DiversityData.csv", sep = ""),
            sep = ";", row.names = F)
