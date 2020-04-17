########################################################################
# NOTES:
# The input data for this file is not part of the GitHub repository,
# but stored elsewhere.
# This file is just used for pre-processing and its final result - the 
# data file fragSAD_predicts_ewers.csv - is provided via GitHub
######################################################################## 

# read long format data file
path2Dropbox <- "C:/Users/May/Dropbox (Privat)/Frag Database (new)/"
infile <- paste(path2Dropbox, "files_datapaper/Long_format_database/fragSAD_long.csv", sep = "")
dat_long <- read.csv(infile, stringsAsFactors = F)
dim(dat_long)
str(dat_long)

# read and adjust files from PREDICTS database

# Mapping fragSAD columns to PREDICTS columns
# fragSAD            PREDICTS
# ------------------------
# frag_id            Site_name
# sample_id          ??? Block maybe
# frag_size_char     Just converted from Habitat_patch_area_square_metres
# frag_size_num      Habitat_patch_area_square_metres --> convert to ha !!!
# sample_eff         Sampling effort --> does rescaling make sense?
# sample_design      Defined by us!
# species            Taxon_name_entered
# abundance          Measurement   

# Caceres 2010 ------------------------------------------------------------
predicts_path <- paste(path2Dropbox, "From PREDICTS/Caceres_et_al_2010.csv", sep = "") 

caceres1 <- read.csv(predicts_path, stringsAsFactors = F)

# Check data set
dim(caceres1)
dim(distinct(caceres1))

# Select and rename columns and filter positive abundances
caceres2 <- caceres1 %>%
   select(Site_name,
          Habitat_patch_area_square_metres,
          Sampling_effort,
          Taxon_name_entered,
          Measurement) %>%
   rename(frag_id       = Site_name,
          frag_size_num = Habitat_patch_area_square_metres,
          sample_eff    = Sampling_effort,
          species       = Taxon_name_entered,
          abundance     = Measurement) %>%
   filter(abundance > 0)
rm(caceres1)

# Add or convert columns
caceres2$dataset_label  <- "Caceres_2010"
caceres2$sample_id      <- 1
caceres2$frag_size_num  <- caceres2$frag_size_num/10000
caceres2$frag_size_char <- as.character(caceres2$frag_size_num)

# Deduce sampling design
# Compare to Table 1 in Caceres et al. 2010
caceres2 %>% 
   select(frag_id, frag_size_num, sample_eff) %>% 
   distinct() %>%
   arrange(frag_size_num)

# different effort in different fragments!!!
caceres2$sample_design <- "pooled"

# sort columns
caceres2a <- caceres2 %>%
   select(names(dat_long))
rm(caceres2)

# sum species abundances in fragments
caceres2b <- caceres2a %>%
   group_by_at(vars(-abundance)) %>%
   summarise(abundance = sum(abundance)) %>% 
   ungroup() %>%
   arrange(frag_id, species)
rm(caceres2a)


# # Ewers 2007 ------------------------------------------------------------
# predicts_path <- path2Dropbox %+% "From PREDICTS/Ewers_et_al_2007.csv" 
# 
# ewers1 <- read.csv(predicts_path, stringsAsFactors = F)
# dim(ewers1)
# dim(distinct(ewers1))
# 
# # No. of species
# length(unique(ewers1$Taxon_name_entered))
# 
# distinct(ewers1,Habitat_as_described)
# 
# ewers2 <- ewers1 %>% 
#    filter(Habitat_as_described != "grassland matrix" & Measurement > 0)
# rm(ewers1)
#    
# length(unique(ewers2$Taxon_name_entered))
# 
# ewers3 <- ewers2 %>%
#    select(Site_name,
#           Habitat_as_described,
#           Habitat_patch_area_square_metres,
#           Sampling_effort,
#           Taxon_name_entered,
#           Measurement) %>%
#    rename(frag_id       = Site_name,
#           frag_size_num = Habitat_patch_area_square_metres,
#           sample_eff    = Sampling_effort,
#           species       = Taxon_name_entered,
#           abundance     = Measurement) %>%
#    filter(abundance > 0)
# 
# length(unique(ewers3$species))
# rm(ewers2)
# 
# # Add or convert columns
# ewers3a <- ewers3 %>%
#    mutate(
#       dataset_label  = "Ewers_2007",
#       sample_id      = 1,
#       frag_size_num  = frag_size_num/10000,
#       frag_size_char = as.character(frag_size_num)
#    )
# 
# # calculate species abundances in each fragment
# ewers3a_abund <- ewers3a %>%
#    group_by(Habitat_as_described, frag_size_char, species) %>%
#    summarise(abundance = sum(abundance))
# 
# # sampling effort in each fragment (sum over sites)
# ewers3a_sample_eff <- ewers3a %>%
#    select(Habitat_as_described, frag_id, frag_size_char, sample_eff) %>%
#    distinct() %>%
#    group_by(Habitat_as_described, frag_size_char) %>%
#    summarise(sample_eff = sum(sample_eff))
# 
# ewers3a_frag <- ewers3a %>%
#    select(Habitat_as_described, frag_size_num, dataset_label, sample_id, frag_size_char) %>%
#    distinct() %>%
#    left_join(ewers3a_sample_eff) %>%
#    arrange(desc(frag_size_num))
# 
# # compare to Table in Ewers et al. 2007
# ewers3a_frag
# outfile <- path2Dropbox %+% "From PREDICTS/Ewers_tab_predicts.csv" 
# write_csv(ewers3a_frag, outfile )
# 
# # Add new fragment IDs
# ewers3a_frag$frag_id <- paste("Fragment" ,1:nrow(ewers3a_frag), sep = "")
# ewers3a_frag$sample_design <- "pooled"
# 
# ewers3b <- left_join(ewers3a_frag, ewers3a_abund)
# 
# length(unique(ewers3$species))
# # just 127 compared to 893 in the paper
# 
# ewers3c <- ewers3b %>%
#    select(names(dat_long))

# Fernandez 2013 ------------------------------------------------------------
predicts_path <- paste(path2Dropbox, "From PREDICTS/Fernandez_&_Simonetti_2013.csv", sep = "") 

fernandez1 <- read.csv(predicts_path, stringsAsFactors = F)

dim(fernandez1)
dim(distinct(fernandez1))

unique(fernandez1$Diversity_metric_type)

fernandez2 <- filter(fernandez1, Diversity_metric_type == "Abundance")
rm(fernandez1)

fernandez2 %>%
   select(Site_name, Habitat_as_described, Habitat_patch_area_square_metres, Sampling_effort) %>%
   distinct()

fernandez3 <- fernandez2 %>%
   select(Site_name,
          Habitat_patch_area_square_metres,
          Sampling_effort,
          Taxon_name_entered,
          Measurement,
          Habitat_as_described) %>%
   rename(frag_id       = Site_name,
          frag_size_num = Habitat_patch_area_square_metres,
          sample_eff    = Sampling_effort,
          species       = Taxon_name_entered,
          abundance     = Measurement) %>%
   filter(abundance > 0 & !is.na( Habitat_as_described))

rm(fernandez2)

# Add or convert columns
fernandez3a <- fernandez3 %>%
   mutate(
      dataset_label  = "Fernandez_2013",
      sample_id      = 1,
      frag_size_num  = frag_size_num/10000,
      frag_size_char = as.character(frag_size_num),
      log_frag_size  = log10(frag_size_num)
   ) %>%
   arrange(frag_size_num, frag_id)

rm(fernandez3)

# Deduce sampling design
fernandez3a %>% 
   select(frag_id, frag_size_num, log_frag_size, sample_eff) %>% distinct

# different effort in different fragments!!!
fernandez3a$sample_design <- "standardized_fragment"

dim(fernandez3a)
dim(distinct(fernandez3a))

# sum species abundances in fragments
fernandez3b <- fernandez3a %>%
   group_by_at(vars(-abundance)) %>%
   summarise(abundance = sum(abundance)) %>% 
   ungroup() %>%
   arrange(frag_id, species)
rm(fernandez3a)

# create sites by species table
fernandez4 <- fernandez3b %>%
   select(frag_id, frag_size_num, species, abundance) %>%
   spread(key = species , value = abundance, fill = 0)

#predicts_path <- path2Dropbox %+% "From PREDICTS/Fernandez_Simonetti_Table_3.csv" 
#write_csv(fernandez4, predicts_path)

fernandez3b_rural <- filter(fernandez3b, Habitat_as_described == "fragments similar in area and habitat characteristics with those of urban area, but surrounded by a rural matrix")
fernandez3b_urban <- filter(fernandez3b, Habitat_as_described == "remnant fragments within an urban matrix")
rm(fernandez3b)

fernandez3b_rural$dataset_label <- "Fernandez_2013_a"
fernandez3b_urban$dataset_label <- "Fernandez_2013_b"

fernandez3c_rural <- fernandez3b_rural %>%
   select(names(dat_long))

fernandez3c_urban <- fernandez3b_urban %>%
   select(names(dat_long))

rm(fernandez3b_rural, fernandez3b_urban)

# Garmendia 2013 ------------------------------------------------------------
predicts_path <- paste(path2Dropbox, "From PREDICTS/Garmendia_et_al_2013.csv", sep = "") 

garmendia1 <- read.csv(predicts_path, stringsAsFactors = F)

# Check data set
dim(garmendia1)
dim(distinct(garmendia1))

garmendia1 %>% distinct(Predominant_land_use)

garmendia1 %>%
   select(Site_name, Predominant_land_use, Habitat_patch_area_square_metres) %>%
   distinct() %>%
   arrange(Predominant_land_use)

unique(garmendia1$Site_name) # less sites than in paper

# Filter for primarey vegetation patches
# Not necessary according to Info from Victor Arroyo-Rdriguez (e-mail from March 22nd 2019)
# garmendia2 <- filter(garmendia1, Predominant_land_use == "Primary vegetation")

# Select and rename columns and filter positive abundances
# filter distinct rows
garmendia3 <- garmendia1 %>%
   select(Site_name,
          Habitat_patch_area_square_metres,
          Sampling_effort,
          Taxon_name_entered,
          Measurement) %>%
   rename(frag_id       = Site_name,
          frag_size_num = Habitat_patch_area_square_metres,
          sample_eff    = Sampling_effort,
          species       = Taxon_name_entered,
          abundance     = Measurement) %>%
   filter(abundance > 0) %>%
   distinct()
rm(garmendia1)

# Add or convert columns
garmendia3a <- garmendia3 %>%
   mutate(
      dataset_label  = "Garmendia_2013",
      sample_id      = 1,
      frag_size_num  = frag_size_num/10000,
      frag_size_char = as.character(frag_size_num)
   ) %>%
   arrange(frag_id, species)
rm(garmendia3)

# Deduce sampling design
garmendia3a %>% select(frag_id, frag_size_num, sample_eff) %>%
   distinct() %>% arrange(frag_size_num)


# set numeric fragment areas to NA 
# and character fragment areas to "continuous" for continuous areas
i_cont <- str_detect(garmendia3a$frag_id, "CF")
garmendia3a$frag_size_num[i_cont] <- NA
garmendia3a$frag_size_char[i_cont] <- "continuous"

# The control sites are in the same continuous forest 
# substitute different sample IDs
garmendia3a$sample_id[i_cont] <- as.integer(str_sub(garmendia3a$frag_id[i_cont],3,4))
garmendia3a$frag_id[i_cont] <- "CF"

# different effort in different fragments!!!
garmendia3a$sample_design <- "standardized_subsamples"

# sort columns
garmendia3b <- garmendia3a %>%
   select(names(dat_long))
rm(garmendia3a)

dim(garmendia3b)
dim(distinct(garmendia3b))

# Stouffer 2011 ------------------------------------------------------------
predicts_path <- paste(path2Dropbox, "From PREDICTS/Stouffer_et_al_2011.csv", sep = "") 
stouffer1 <- read.csv(predicts_path, stringsAsFactors = F)
dim(stouffer1)
dim(distinct(stouffer1))

stouffer1 %>% distinct(Predominant_land_use)
stouffer1 %>% distinct(Habitat_patch_area_square_metres)

stouffer2 <- stouffer1 %>%
   select(Site_name,
          Habitat_patch_area_square_metres,
          Sampling_effort,
          Taxon_name_entered,
          Measurement) %>%
   rename(frag_id       = Site_name,
          frag_size_num = Habitat_patch_area_square_metres,
          sample_eff    = Sampling_effort,
          species       = Taxon_name_entered,
          abundance     = Measurement) %>%
   filter(abundance > 0)
rm(stouffer1)

# Add or convert columns
stouffer2a <- stouffer2 %>%
   mutate(
      dataset_label  = "Stouffer_2011",
      sample_id      = 1,
      frag_size_num  = frag_size_num/10000,
      frag_size_char = as.character(frag_size_num)
   )
rm(stouffer2)

# Deduce sampling design
stouffer2a %>% select(frag_id, frag_size_num, sample_eff) %>%
   distinct() %>% arrange(frag_size_num)

# different effort in different fragments!!!
stouffer2a$sample_design <- "pooled"

stouffer2b <- stouffer2a %>%
   select(names(dat_long))
stouffer2b$frag_id <- as.character(stouffer2b$frag_id)
rm(stouffer2a)

# Combine and save studies ----------------------------------
# exclude ewers for the moment
predicts_dat <- bind_rows(caceres2b,
                          fernandez3c_rural, fernandez3c_urban,
                          garmendia3b,
                          stouffer2b)

# Remember to standardize sampling effort within studies!!! 
predicts_dat2 <- predicts_dat %>% 
   group_by(dataset_label) %>%
   mutate(sample_eff = sample_eff/min(sample_eff))

predicts_dat2 %>% 
   summarise(min_sample_eff = min(sample_eff))

# View(predict_dat2 %>% 
#    group_by(dataset_label) %>%
#    select(dataset_label, frag_id, sample_eff) %>% distinct())

fragsad_predicts <- bind_rows(dat_long, predicts_dat)
sum(duplicated(fragsad_predicts))

# read long format data file
outfile <- paste(path2Dropbox, "files_datapaper/Long_format_database/fragSAD_and_predicts.csv", sep = "")
write_csv(fragsad_predicts, outfile)
