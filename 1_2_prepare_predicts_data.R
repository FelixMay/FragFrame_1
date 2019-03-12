# read long format data file
infile <- path2Dropbox %+% "files_datapaper/Long_format_database/fragSAD_long.csv"
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
predicts_path <- path2Dropbox %+% "From PREDICTS/Caceres_et_al_2010.csv" 
caceres1 <- read.csv(predicts_path, stringsAsFactors = F)

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

# Add or convert columns
caceres2$dataset_label  <- "Caceres_2010"
caceres2$sample_id      <- 1
caceres2$frag_size_num  <- caceres2$frag_size_num/10000
caceres2$frag_size_char <- as.character(caceres2$frag_size_num)

# Deduce sampling design
caceres2 %>% select(frag_id, frag_size_num, sample_eff) %>% distinct

# different effort in different fragments!!!
caceres2$sample_design <- "pooled"

caceres2 <- caceres2 %>%
   select(names(dat_long)) %>%
   mutate(sample_eff = sample_eff/min(sample_eff))

# Ewers 2007 ------------------------------------------------------------
predicts_path <- path2Dropbox %+% "From PREDICTS/Ewers_et_al_2007.csv" 
ewers1 <- read.csv(predicts_path, stringsAsFactors = F)

distinct(ewers1,Habitat_as_described)
ewers1 <- ewers1 %>% filter(Habitat_as_described != "grassland matrix")

ewers2 <- ewers1 %>%
   select(Site_name,
          Habitat_patch_area_square_metres,
          Sampling_effort,
          Taxon,
          Measurement) %>%
   rename(frag_id       = Site_name,
          frag_size_num = Habitat_patch_area_square_metres,
          sample_eff    = Sampling_effort,
          species       = Taxon,
          abundance     = Measurement) %>%
   filter(abundance > 0)

# Add or convert columns
ewers2 <- ewers2 %>%
   mutate(
      dataset_label  = "Ewers_2007",
      sample_id      = 1,
      frag_size_num  = frag_size_num/10000,
      frag_size_char = as.character(frag_size_num)
   )

# calculate species abundances in each fragment
ewers_abund <- ewers2 %>%
   group_by(frag_size_char, species) %>%
   summarise(abundance = sum(abundance))

ewers_sample_eff <- ewers2 %>%
   select(frag_id, frag_size_char, sample_eff) %>%
   distinct() %>%
   group_by(frag_size_char) %>%
   summarise(sample_eff = sum(sample_eff))

ewers_frag <- ewers2 %>%
   select(frag_size_num, dataset_label, sample_id, frag_size_char) %>%
   distinct() %>%
   left_join(ewers_sample_eff) %>%
   arrange(desc(frag_size_num))
ewers_frag$frag_id <- paste("Site",1:nrow(ewers_frag), sep = "")
ewers_frag$sample_design <- "pooled"

ewers3 <- left_join(ewers_frag, ewers_abund)

ewers3 <- ewers3 %>%
   select(names(dat_long))

# Fernandez 2013 ------------------------------------------------------------
predicts_path <- path2Dropbox %+% "From PREDICTS/Fernandez_&_Simonetti_2013.csv" 
fernandez1 <- read.csv(predicts_path, stringsAsFactors = F)

distinct(fernandez1, Habitat_as_described)

fernandez2 <- fernandez1 %>%
   select(Site_name,
          Habitat_patch_area_square_metres,
          Sampling_effort,
          Taxon,
          Measurement,
          Habitat_as_described) %>%
   rename(frag_id       = Site_name,
          frag_size_num = Habitat_patch_area_square_metres,
          sample_eff    = Sampling_effort,
          species       = Taxon,
          abundance     = Measurement) %>%
   filter(abundance > 0)

# Add or convert columns
fernandez2 <- fernandez2 %>%
   mutate(
      dataset_label  = "Fernandez_2013",
      sample_id      = 1,
      frag_size_num  = frag_size_num/10000,
      frag_size_char = as.character(frag_size_num)
   )

# Deduce sampling design
fernandez2 %>% select(frag_id, frag_size_num, sample_eff) %>% distinct

# different effort in different fragments!!!
fernandez2$sample_design <- "pooled"


fernandez2_rural <- filter(fernandez2, Habitat_as_described == "fragments similar in area and habitat characteristics with those of urban area, but surrounded by a rural matrix")
fernandez2_urban <- filter(fernandez2, Habitat_as_described == "remnant fragments within an urban matrix")

fernandez2_rural$dataset_label <- "Fernandez_2013_a"
fernandez2_urban$dataset_label <- "Fernandez_2013_b"

fernandez2_rural <- fernandez2_rural %>%
   select(names(dat_long))

fernandez2_urban <- fernandez2_urban %>%
   select(names(dat_long))


# Garmendia 2013 ------------------------------------------------------------
predicts_path <- path2Dropbox %+% "From PREDICTS/Garmendia_et_al_2013.csv" 
garmendia1 <- read.csv(predicts_path, stringsAsFactors = F)

garmendia1 %>% distinct(Predominant_land_use)

garmendia1 <- filter(garmendia1, Predominant_land_use == "Primary vegetation")

garmendia2 <- garmendia1 %>%
   select(Site_name,
          Habitat_patch_area_square_metres,
          Sampling_effort,
          Taxon,
          Measurement) %>%
   rename(frag_id       = Site_name,
          frag_size_num = Habitat_patch_area_square_metres,
          sample_eff    = Sampling_effort,
          species       = Taxon,
          abundance     = Measurement) %>%
   filter(abundance > 0)

# Add or convert columns
garmendia2 <- garmendia2 %>%
   mutate(
      dataset_label  = "Garmendia_2013",
      sample_id      = 1,
      frag_size_num  = frag_size_num/10000,
      frag_size_char = as.character(frag_size_num)
   )

# Deduce sampling design
garmendia2 %>% select(frag_id, frag_size_num, sample_eff) %>%
   distinct() %>% arrange(frag_size_num)

# different effort in different fragments!!!
garmendia2$sample_design <- "standardized"

garmendia2 <- garmendia2 %>%
   select(names(dat_long))

# Stouffer 2011 ------------------------------------------------------------
predicts_path <- path2Dropbox %+% "From PREDICTS/Stouffer_et_al_2011.csv" 
stouffer1 <- read.csv(predicts_path, stringsAsFactors = F)

stouffer1 %>% distinct(Predominant_land_use)

stouffer2 <- stouffer1 %>%
   select(Site_name,
          Habitat_patch_area_square_metres,
          Sampling_effort,
          Taxon,
          Measurement) %>%
   rename(frag_id       = Site_name,
          frag_size_num = Habitat_patch_area_square_metres,
          sample_eff    = Sampling_effort,
          species       = Taxon,
          abundance     = Measurement) %>%
   filter(abundance > 0)

# Add or convert columns
stouffer2 <- stouffer2 %>%
   mutate(
      dataset_label  = "Stouffer_2011",
      sample_id      = 1,
      frag_size_num  = frag_size_num/10000,
      frag_size_char = as.character(frag_size_num)
   )

# Deduce sampling design
stouffer2 %>% select(frag_id, frag_size_num, sample_eff) %>%
   distinct() %>% arrange(frag_size_num)

# different effort in different fragments!!!
stouffer2$sample_design <- "pooled"

stouffer2 <- stouffer2 %>%
   select(names(dat_long))
stouffer2$frag_id <- as.character(stouffer2$frag_id)

# Combine and save studies
predicts_dat <- bind_rows(caceres2, ewers2, fernandez2_rural, fernandez2_urban, garmendia2, stouffer2)

