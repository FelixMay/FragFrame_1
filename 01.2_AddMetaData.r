setwd(path2Dropbox %+% "good_datasets/")


# ---------------------------------------------------
### join meta data extracted from files with manual extensions
meta_df <- read.xlsx(path2Dropbox %+% "_Mario data curating/clean metadata matrix_version_2017_December_18.xlsx", sheetIndex = 2, startRow = 1,  header = T)

# ---------------------------------------------------
### group taxa levels
meta_df$taxa <- as.character(meta_df$taxa)
meta_df$taxa[meta_df$taxa %in% c("amphibians","amphibians, reptiles", "reptiles", "reptiles, molluscs")] <- "amphibians & reptiles"
meta_df$taxa[meta_df$taxa %in% c("birds")] <- "birds"
meta_df$taxa[meta_df$taxa %in% c("mammals")] <- "mammals"
meta_df$taxa[meta_df$taxa %in% c("plants")] <- "plants"
meta_df$taxa[meta_df$taxa %in% c("arachnids", "insects", "arthropods")] <- "invertebrates"

meta_df$taxa <- factor(meta_df$taxa) 

# ---------------------------------------------------
### make levels consistent
## relevel matrix category
meta_df$matrix.filter.category <- as.character(meta_df$matrix.filter.category)
meta_df$matrix.filter.category[meta_df$matrix.filter.category == "intermediate"] <- "medium filter"
meta_df$matrix.filter.category <- factor(meta_df$matrix.filter.category, levels=c("light filter", "medium filter", "harsh filter")) 

## relevel time.since.fragmentation
meta_df$time.since.fragmentation <- as.character(meta_df$time.since.fragmentation)
meta_df$time.since.fragmentation[meta_df$time.since.fragmentation == "Intermediate (20-100 years)"] <- "intermediate (20-100 years)"
meta_df$time.since.fragmentation[meta_df$time.since.fragmentation == "Recent (less than 20 years)"] <- "recent (less than 20 years)"
meta_df$time.since.fragmentation <- factor(meta_df$time.since.fragmentation, levels=c("recent (less than 20 years)", "intermediate (20-100 years)", "long (100+ years)")) 

# ---------------------------------------------------
### save table
write.csv(meta_df, file = paste(path2temp, "metaData.csv", sep = ""))

##############################
### RESTERAMPE
# ---------------------------------------------------
### extract meta data from xls-sheets
# filenames <- list.files(pattern="*.xls*", full.names = F)
# filenames2 <- sapply(strsplit(filenames, split = "[.]"), "[[", 1)
# 
# meta_list <- list()
# 
# for (i in 1:length(filenames)){
#    
#    print(filenames[i])
#    
#    df <- try(read.xlsx(filenames[i], sheetIndex = 2, rowIndex=1:23, colIndex=1:2,  header = F)) # specify row and column indices to avoid reading empty cols and rows which would produce errors
#    if (!inherits(df, "try-error")){
#       if(ncol(df)>1){
#          dat_meta <- as.data.frame(matrix(data=NA,ncol=nrow(df),nrow=1))
#          dat_meta[1,] <- t(df[,2])
#          names(dat_meta) <- df[,1]
#          
#          meta_list[[filenames2[i]]] <- data.frame(Case.ID = filenames2[i],
#                                                   taxa=dat_meta$taxa.class,
#                                                   country=dat_meta$location.country,
#                                                   continent=dat_meta$location.continent,
#                                                   biome=dat_meta$location.biome,
#                                                   fragment.biome = dat_meta$fragment.biome,
#                                                   matrix.biome=dat_meta$matrix.biome,
#                                                   fragment.veg=dat_meta$fragment.vegetation,
#                                                   matrix.veg=dat_meta$matrix.vegetation,
#                                                   sampling.effort = dat_meta$sampling.effort.measure)
#          meta_list[[filenames2[i]]] <- sapply(meta_list[[filenames2[i]]],as.character)
#       }
#    }
# }
# 
# meta_df <- as.data.frame(do.call(rbind, meta_list) )
