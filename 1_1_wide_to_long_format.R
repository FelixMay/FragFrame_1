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
   if (length(duplicated_species) > 0)
      warning(paste("Duplicated species", duplicated_species, sep = ", ")," in ", filename)
   
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
   
   return(dat1)
}    

################################################################################
# execution of script

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

