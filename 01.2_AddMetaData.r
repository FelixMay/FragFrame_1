setwd(path2Dropbox %+% "good_datasets/")

filenames <- list.files(pattern="*.xls*", full.names = F)
filenames2 <- sapply(strsplit(filenames, split = "[.]"), "[[", 1)

meta_list <- list()

for (i in 1:length(filenames)){
   
   print(filenames[i])
   
   df <- try(read.xlsx(filenames[i], sheetIndex = 2, startRow = 1,  header = F))
   if (!inherits(df, "try-error")){
      if(ncol(df)>1){
         dat_meta <- as.data.frame(matrix(data=NA,ncol=nrow(df),nrow=1))
         dat_meta[1,] <- t(df[,-1])
         names(dat_meta) <- df[,1]
         
         meta_list[[filenames2[i]]] <- data.frame(Study.ID = filenames2[i],
                                                  taxa=dat_meta$taxa.class,
                                                  country=dat_meta$location.country,
                                                  continent=dat_meta$location.continent,
                                                  biome=dat_meta$location.biome)
         meta_list[[filenames2[i]]] <- sapply(meta_list[[filenames2[i]]],as.character)
      }
   }
}

meta_df <- as.data.frame(do.call(rbind, meta_list) )

write.csv(meta_df, file = paste(path2temp, "metaData.csv", sep = ""))
