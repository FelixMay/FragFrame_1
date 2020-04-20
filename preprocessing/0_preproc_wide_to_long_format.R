########################################################################
# NOTES:
# The input data for this file is not part of the GitHub repository,
# but stored elsewhere.
# This file is just used for pre-processing and its final result - the 
# data file fragSAD_predicts_ewers.csv - is provided via GitHub
######################################################################## 

# Read files, transpose and switch to long format

read_data_files <- function(filename){
   
   print(filename)

   path2infile <- paste(path2Dropbox, filename, sep = "")
   
   # 1. Abundance data
   dat_abund <- read.table(path2infile, sep = ",", row.names = NULL,
                           skip = 4, stringsAsFactors = F)
   
   ## delete columns with only NAs
   na_per_col <- apply(dat_abund, 2, function(x) sum(is.na(x)))
   dat_abund <- dat_abund[, na_per_col < nrow(dat_abund)]
   
   n_col <- ncol(dat_abund) - 1 # No. of columns with abundance data
   
   ## check duplicates in first column
   duplicated_species <- dat_abund[[1]][duplicated(dat_abund[[1]])]
   if (length(duplicated_species) > 0){
      names_duplicates <- paste(duplicated_species, collapse = ", ")
      print(paste("Duplicated species in", filename, names_duplicates, sep = " "))
   }
   
   spec_names <- as.character(dat_abund[[1]])
   dat_abund <- dat_abund[,-1]
   dat_abund_t <- as.data.frame(t(dat_abund))
   names(dat_abund_t) <- spec_names
   
   ## sort by species name for better spotting of duplicates
   dat_abund_t <- dat_abund_t[, order(names(dat_abund_t))]
   
   ## Check for non-numeric abundances
   is_num <- sapply(dat_abund_t, is.numeric)
   if (!all(is_num))
      print(paste("Non-numeric abundances in", filename, sep = " "))
   
   # 2. Sampling effort
   sample_eff <- read.table(path2infile, sep = ",", row.names = 1, nrows = 1,
                            stringsAsFactors = F)
   
   ## remove non-numeric entries
   sample_eff <- as.numeric(sample_eff[,1:n_col])
   
   if (any(is.na(sample_eff)))
      print(paste("NA values in sampling effort in", filename, sep = " "))
   
   if (max(sample_eff) > 1)
      print(paste("Unequal sampling effort in", filename))
   
   # 3. Fragment ID
   frag_id <- read.table(path2infile, sep = ",", row.names = 1, nrows = 1,
                         skip = 1, stringsAsFactors = F)
   frag_id <- as.character(frag_id[, 1:n_col])
   
   # 4. sample ID
   sample_id <- read.table(path2infile, sep = ",", row.names = 1, nrows = 1,
                                skip = 2, stringsAsFactors = F)
   sample_id <- as.character(sample_id[, 1:n_col])
   
   # 5. Fragment size
   frag_size <- read.table(path2infile, sep = ",", row.names = 1, nrows = 1,
                             skip = 3, stringsAsFactors = F)
   frag_size_char <- as.character(frag_size[, 1:n_col])
   frag_size_num <- as.numeric(frag_size_char)
   
   if (any(is.na(frag_size_num)))
      print(paste("NA values in fragment sizes in", filename, sep = " "))
   
   # combine fragment-data
   dat1 <- data.frame(frag_id,
                      sample_id,
                      frag_size_char,
                      frag_size_num,
                      sample_eff,
                      stringsAsFactors = F)
   
   # determine sampling design
   sample_eff_per_frag <- tapply(dat1$sample_eff, dat1$frag_id, sum)
   plot_id <- paste(dat1$frag_id, dat1$sample_id, sep = "_")
   sample_eff_per_plot <- tapply(dat1$sample_eff, plot_id, sum)
   range_sample_eff_frag <- max(sample_eff_per_frag) - min(sample_eff_per_frag)   
   range_sample_eff_plot <- max(sample_eff_per_plot) - min(sample_eff_per_plot) 
   
   n_subsamples <- table(dat1$frag_id)
   
   if (max(n_subsamples) == 1 && range_sample_eff_frag == 0 && range_sample_eff_plot == 0){
      dat1$sample_design <- "standardized_fragment"
   } else {
      if (range_sample_eff_plot == 0) {
         dat1$sample_design <- "standardized_subsamples"
      } else {
         dat1$sample_design <- "pooled"
      }
   }
   
   # combine fragment and abundance data
   dat1 <- cbind(dat1, dat_abund_t)
   
   # prepare output date
   outfile <- str_split(filename, pattern = "_\\(")[[1]][1]
   path2outfile <- paste(path2Dropbox, "Sites_by_species_format/", outfile, ".csv", sep = "")
   write_csv(dat1, path2outfile)
   
   # convert to long format: one column with species name and one with species abundance
   dat2 <- dat1 %>%
      gather(-(frag_id:sample_design), key = "species", value = "abundance") %>%
      filter(abundance > 0) %>%
      arrange(frag_id, sample_id, species)
   
   # labels <- str_split(outfile, "_")
   # dat2$dataset_id <- labels[[1]][1]
   dat2$dataset_label <- outfile
   
   dat2 <- select(dat2, dataset_label, everything())
   
   # check for non-integer abundances
   if (!is.integer(dat2$abundance))
      print(paste("Non-integer abundances in", filename, sep = " "))
   
   cat("\n")
   
   return(dat2)
}    

################################################################################
# execution of script
path2Dropbox <- "C:/Users/May/Dropbox (Privat)/Frag Database (new)/files_datapaper/"

filenames <- list.files(path =  path2Dropbox, pattern="*.csv", full.names = F)

dat_list <- lapply(filenames, read_data_files)

dat_long <- bind_rows(dat_list)

summary(dat_long)

#dat_long %>% filter(species == "entity.plot")
dat_long %>% filter(frag_size_char == "Continuous")

# prepare output date
path2outfile <- paste(path2Dropbox, "Long_format_database/fragSAD_long.csv", sep = "")
write_csv(dat_long, path2outfile)



