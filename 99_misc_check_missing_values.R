infile <- path2Dropbox %+% "/files_datapaper/Analysis/2_biodiv_frag_fcont_10_mabund_as_is.csv"
frag_div <- read.csv(infile, stringsAsFactors = F)

summary(frag_div)

na_S_std_1 <- frag_div %>%
   filter(is.na(S_std_1)) 
summary(na_S_std_1)
# seems to be OK --> N_std is always smaller than 5 and than rarefaction does not work

na_S_std_2 <- frag_div %>%
   filter(is.na(S_std_2)) 
summary(na_S_std_2)
# same just for the pooled sampling designs (in the other ones rarefaction is not used)

na_S_PIE <- frag_div %>%
   filter(is.na(S_PIE)) 
summary(na_S_PIE)
dim(na_S_PIE)
select(na_S_PIE, N, S_obs)
# S == N!

na_S_n <- frag_div %>%
   filter(is.na(S_n)) 
na_S_n
summary(na_S_n)

frag_div %>% filter(dataset_label == "Ewers_2007")

hist(unique(frag_div$cov_base))

# check base coverage
frag_div_cov_base <- frag_div %>%
   select(dataset_label, cov_base) %>%
   distinct() %>%
   arrange(cov_base) 

frag_div_cov <- frag_div %>%
   select(dataset_label, coverage, cov_base) %>%
   distinct() %>%
   arrange(cov_base) 
plot(coverage ~ cov_base, data = frag_div_cov)
hist(frag_div_cov$coverage)
