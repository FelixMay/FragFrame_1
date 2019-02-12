# https://stackoverflow.com/questions/23474729/convert-object-of-class-dist-into-data-frame-in-r
dist_to_dataframe <- function(inDist) {
   if (class(inDist) != "dist") stop("wrong input type")
   A <- attr(inDist, "Size")
   B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
   if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
   if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
   data.frame(
      row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
      col = rep(B[-length(B)], (length(B)-1):1),
      value = as.vector(inDist),
      stringsAsFactors = F
   )
}

############################################
# Function to calculate beta-diversity partitioning for every
# pair of fragments

get_beta_part <- function(data_set){
   
   print(data_set$dataset_label[1])
   
   # sum abundances in same fragments
   dat_abund <- data_set %>%
      group_by(frag_id, species) %>%
      summarise(abundance = sum(abundance))
   
   dat_wide <- dat_abund %>% spread(key = species, value = abundance, fill = 0)
   
   # add rownames for betadivpart
   class(dat_wide) <- "data.frame"
   rownames(dat_wide) <- dat_wide$frag_id
   
   # # prepare output data
   # dat_sample <- data_set %>%
   #    select(-species, -abundance) %>%
   #    distinct()
   # 
   # dat_sample_eff <- dat_sample %>%
   #    group_by(frag_id) %>%
   #    summarise(sample_eff = sum(sample_eff)) 
   # 
   # dat_frag <- dat_sample %>% 
   #    select(-sample_id, -sample_eff) %>%
   #    distinct() %>%
   #    left_join(dat_sample_eff)
   
   # dat_wide <- dat_frag %>%
   #    select(frag_id, sample_eff) %>%
   #    left_join(dat_wide)
   
   #############################################################################
   ### Beta-diversity partitioning
   ### Estimate species turnover that is due to replacement (in contrast to nestedness)
   
   J_qF <- beta.div.comp(dat_wide[,-1], coef = "J", quant = F)
   S_qF <- beta.div.comp(dat_wide[,-1], coef = "S", quant = F)
   BJ_qF <- beta.div.comp(dat_wide[,-1], coef = "BJ", quant = F)
   BS_qF <- beta.div.comp(dat_wide[,-1], coef = "BS", quant = F)
   
   J_qT <- beta.div.comp(dat_wide[,-1], coef = "J", quant = T)
   S_qT <- beta.div.comp(dat_wide[,-1], coef = "S", quant = T)
   BJ_qT <- beta.div.comp(dat_wide[,-1], coef = "BJ", quant = T)
   BS_qT <- beta.div.comp(dat_wide[,-1], coef = "BS", quant = T)
   
   beta_div_tab <- dist_to_dataframe(J_qF$repl)
   names(beta_div_tab)[1:3] <- c("frag_x", "frag_y", "J_qF_repl")
   
   beta_div_tab$J_qF_rich <- J_qF$rich
   
   beta_div_tab$S_qF_repl <- S_qF$repl
   beta_div_tab$S_qF_rich <- S_qF$rich
   
   beta_div_tab$BJ_qF_repl <- BJ_qF$repl
   beta_div_tab$BJ_qF_rich <- BJ_qF$rich
   
   beta_div_tab$BS_qF_repl <- BS_qF$repl
   beta_div_tab$BS_qF_rich <- BS_qF$rich
   
   beta_div_tab$J_qT_repl <- J_qF$repl
   beta_div_tab$J_qT_rich <- J_qF$rich
   
   beta_div_tab$S_qT_repl <- S_qT$repl
   beta_div_tab$S_qT_rich <- S_qT$rich
   
   beta_div_tab$BJ_qT_repl <- BJ_qT$repl
   beta_div_tab$BJ_qT_rich <- BJ_qT$rich
   
   beta_div_tab$BS_qT_repl <- BS_qT$repl
   beta_div_tab$BS_qT_rich <- BS_qT$rich
   
   # add fragment areas
   dat_frag <- data_set %>%
      select(frag_id, frag_size_char, frag_size_num) %>%
      distinct()

   # dat_frag$site_label <- paste("Site", 1:nrow(dat_frag), sep = "")
   beta_div_tab <- beta_div_tab %>%
      left_join(dat_frag, by = c("frag_x" = "frag_id")) %>%
      left_join(dat_frag, by = c("frag_y" = "frag_id")) %>%
      mutate(diff_area = frag_size_num.y - frag_size_num.x,
             log10_ratio_area = log10(frag_size_num.y/frag_size_num.x)
             )
   
   beta_div_tab$dataset_id <- data_set$dataset_id[1]
   beta_div_tab$dataset_label <- data_set$dataset_label[1]
   beta_div_tab$sample_design <- data_set$sample_design[1]
   
   beta_div_tab <- beta_div_tab %>% 
      select(dataset_id,
             dataset_label,
             sample_design,
             frag_x,
             frag_y,
             frag_size_char.x,
             frag_size_num.x,
             frag_size_char.y,
             frag_size_num.y,
             diff_area,
             log10_ratio_area,
             everything())

   return(beta_div_tab)
}    

################################################################################
# execution of script

# read long format data file

infile <- path2Dropbox %+% "files_datapaper/Long_format_database/fragSAD_long.csv"
dat_long1 <- read.csv(infile, stringsAsFactors = F)
dat_long2 <- read_csv(infile, quoted_na = FALSE,
                      col_types = list(col_integer(),
                                       col_character(),
                                       col_character(),
                                       col_integer(),
                                       col_character(),
                                       col_double(),
                                       col_double(),
                                       col_character(),
                                       col_character(),
                                       col_double()))
                                       
                        
dim(dat_long1)
dim(dat_long2)
# str(dat_long)

# head(dat_long)

# data_set <- dat_long %>% filter(dataset_id == "111")

dat_long1 %>% filter(dataset_id == "111" & frag_id == "NA")
dim(dat_long1 %>% filter(dataset_id == "111" & is.na(frag_id)))

?dat_long2 %>% filter(dataset_id == "111" & frag_id == "NA")
dat_long2 %>% filter(dataset_id == "111" & is.na(frag_id))
is.na(data_set$frag_id)

# purrr version
out1 <- dat_long %>%
   split(.$dataset_id) %>%
   map_dfr(get_beta_part)




