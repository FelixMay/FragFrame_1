########################################################################
# NOTES:
# The input data for this file is not part of the GitHub repository,
# but stored elsewhere.
# This file is just used for pre-processing and its final result - the 
# data file fragSAD_predicts_ewers.csv - is provided via GitHub
######################################################################## 

path2Dropbox <- "C:/Users/May/Dropbox (Privat)/Frag Database (new)/"

infile <- paste(path2Dropbox, "From PREDICTS/Ewers_et_al_2007/HRFFPbeetledata_site.csv", sep = "")
site <- read.csv(infile, stringsAsFactors = F)
names(site)

site <- site %>%
   mutate(Area = 10^log10Area - 1) %>%
   arrange(Site)

n_frag <- length(unique(site$Area))

# Add unique fragment IDs
frag <- site %>% 
   select(Area) %>%
   distinct(Area) %>%
   arrange(Area) %>%
   mutate(frag_id = paste("frag", 1:n_frag, sep = ""))
frag$frag_id[frag$Area > 1000000] <- "continuous1"

site <- left_join(site, frag)
          
sort(unique(site$Area))

# Add sampling effort
infile <- paste(path2Dropbox, "From PREDICTS/Ewers_et_al_2007/HRFFPbeetledata_sample_days.csv", sep = "")
sample_days <- read.csv(infile, stringsAsFactors = F)
sample_days <- arrange(sample_days, NewSite)
head(sample_days)
tail(sample_days)
dim(sample_days)

str(sample_days)
any(duplicated(sample_days))

site$sample_eff <- sample_days$SampDays[1:nrow(site)]
site$Site2      <- sample_days$NewSite[1:nrow(site)]
site <- select(site, Site, Site2, everything())
head(site)
tail(site)

# read abundance data
infile <- paste(path2Dropbox, "From PREDICTS/Ewers_et_al_2007/HRFFPbeetledata_rawdata.csv", sep = "")
raw_data <- read.csv(infile, row.names = 1, stringsAsFactors = F)
dim(raw_data)

raw_data_t <- as.data.frame(t(raw_data))
raw_data_t <- rownames_to_column(raw_data_t,"Site")
dim(raw_data_t)

# switch to long format
raw_data_long <- raw_data_t %>%
   gather(key = species, value = abundance, -Site) %>%
   filter(abundance > 0) %>% arrange(Site, desc(abundance))

site_abund <- left_join(site, raw_data_long)

site_abund %>% 
   select(Site, sample_eff) %>%
   distinct()

# filter samples within forest fragments
site_abund2 <- site_abund %>%
   filter(DistCode < 0)

# adjust to structure of exiting data base
infile <- paste(path2Dropbox, "files_datapaper/Long_format_database/fragSAD_and_predicts.csv", sep = "")
dat_long <- read.csv(infile, stringsAsFactors = F)
dim(dat_long)
str(dat_long)
dat_long$sample_id <- as.character(dat_long$sample_id)

site_abund3 <- site_abund2 %>%
   rename(sample_id = Site,
          frag_size_num = Area) %>%
   mutate(dataset_label = "Ewers_2007",
          frag_size_char = as.character(frag_size_num),
          sample_design  = "pooled") %>%
   select(names(dat_long))
head(site_abund3)

# Remember to standardize sampling effort within studies!!!
site_abund3 <- site_abund3 %>%
   mutate(sample_eff = sample_eff/min(sample_eff) )

fragsad_predicts_ewers <- bind_rows(dat_long, site_abund3)

# read long format data file
outfile <- paste(path2Dropbox, "files_datapaper/Long_format_database/fragSAD_predicts_ewers.csv", sep = "")
write_csv(fragsad_predicts_ewers, outfile)

sum(duplicated(fragsad_predicts_ewers))
