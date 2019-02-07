# Read files, transpose and switch to long format

read_data_files <- function(filename){
   
   print(filename)

   path2infile <- path2Dropbox %+% "files_datapaper/" %+% filename
   
   # Sampling effort
   sample_eff <- read.table(path2infile, sep = ",", row.names = 1, nrows = 1,
                            stringsAsFactors = F)
   
   # remove non-numeric entries
   is_num <- sapply(sample_eff, is.numeric)
   sample_eff <- as.numeric(sample_eff[,is_num])
   
   n_col <- length(sample_eff)
   
   if (max(sample_eff) > 1) print(paste("Unequal sampling effort in", filename))
   
   # entity_id_orig
   entity_id_orig <- read.table(path2infile, sep = ",", row.names = 1, nrows = 1,
                                skip = 1, stringsAsFactors = F)
   entity_id_orig <- as.character(entity_id_orig[1, 1:n_col])
   
   # entity_id_plot
   entity_id_plot <- read.table(path2infile, sep = ",", row.names = 1, nrows = 1,
                                skip = 2, stringsAsFactors = F)
   entity_id_plot <- as.numeric(entity_id_plot[1, 1:n_col])
   
   # entity_size
   entity_size <- read.table(path2infile, sep = ",", row.names = 1, nrows = 1,
                             skip = 3, stringsAsFactors = F)
   entity_size <- as.numeric(entity_size[1, 1:n_col])
   
   # abundance data
   dat_abund <- read.table(path2infile, sep = ",", row.names = NULL,
                           skip = 4, stringsAsFactors = F)
   
   # check duplicates in first column
   duplicated_species <- dat_abund[[1]][duplicated(dat_abund[[1]])]
   if (length(duplicated_species) > 0){
      names_duplicates <- paste(duplicated_species, collapse = ", ")
      print(paste("Duplicated species in", filename, names_duplicates, sep = " "))
   }
     
   spec_names <- as.character(dat_abund[[1]])
   dat_abund <- dat_abund[, 2:(n_col + 1)]
   dat_abund_t <- as.data.frame(t(dat_abund))
   names(dat_abund_t) <- spec_names
   
   # sort by species name for better spotting of duplicates
   dat_abund_t <- dat_abund_t[, order(names(dat_abund_t))]
   
   # combine fragment-data
   dat1 <- data.frame(entity_id_orig, entity_id_plot, entity_size, sample_eff)
   
   # determine sampling design
   sample_eff_per_frag <- tapply(dat1$sample_eff, dat1$entity_id_orig, sum)
   plot_id <- paste(dat1$entity_id_orig, dat1$entity_id_plot, sep = "_")
   sample_eff_per_plot <- tapply(dat1$sample_eff, plot_id, sum)
   range_sample_eff_frag <- max(sample_eff_per_frag) - min(sample_eff_per_frag)   
   range_sample_eff_plot <- max(sample_eff_per_plot) - min(sample_eff_per_plot)  
   
   if (range_sample_eff_frag == 0 && range_sample_eff_plot == 0){
      dat1$sample_design <- "standardized"
   } else {
      if (range_sample_eff_plot == 0) {
         dat1$sample_design <- "subsamples_in_frag"
      } else {
         dat1$sample_design <- "pooled"
      }
   }
   
   # combine fragment and abundance data
   dat1 <- cbind(dat1, dat_abund_t)
   
   # prepare output date
   outfile <- str_split(filename, pattern = "\\s\\(")[[1]][1]
   path2outfile <- path2Dropbox %+% "files_datapaper/Sites_by_species_format/" %+% outfile %+% ".csv"
   write_csv(dat1, path2outfile)
   
   cat("\n")
   
   # convert to long format: one column with species name and one with species abundance
   dat2 <- dat1 %>%
      gather(-(entity_id_orig:sample_design), key = "species", value = "abundance") %>%
      filter(abundance > 0) %>%
      arrange(entity_id_orig, entity_id_plot, species)
   
   labels <- str_split(outfile, "_")
   dat2$dataset_id <- labels[[1]][1]
   dat2$dataset_label <- labels[[1]][2]
   
   dat2 <- select(dat2, dataset_id, dataset_label, everything())
   
   return(dat2)
}    

################################################################################
# execution of script

filenames <- list.files(path =  path2Dropbox %+% "files_datapaper/",
                        pattern="*.csv", full.names = F)

# Files with errors
problem_files <- c("17_Bragagnolo et al. 2007 (harvestmen_in_Brazil).csv",
                   "50_Gavish et al. 2012 (spiders-Galon_in_Israel).csv",
                   "Kapoor 2008 (spiders_in_India).csv")

good_files <- filenames[!filenames %in% problem_files]

length(good_files)
dat_list <- lapply(good_files[80:99], read_data_files)

# #dat_list <- list()
# 
# for (i in 1:length(filenames)){
#    temp <- try(read_data_files(filenames[i]))
#    if (!inherits(temp, "try-error")){
#       dat_list[[filenames[i]]] <- temp
#    }
# }

dat_long <- bind_rows(dat_list)

summary(dat_long)

dat_long %>% filter(species == "entity.size.rank")
unique(dat_long$species)

write.table(div_df_nomatrix, file = paste(path2temp, "DiversityData.csv", sep = ""),
            sep = ",", row.names = F)

