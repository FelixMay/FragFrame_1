infile <- path2Dropbox %+% "/files_datapaper/Analysis/2_biodiv_frag_fcont_10_mabund_as_is.csv"
frag_div <- read.csv(infile, stringsAsFactors = F)

# Check the abundance distribution in fragments where S_PIE > S_obs
infile <- path2Dropbox %+% "files_datapaper/Long_format_database/fragSAD_predicts_ewers.csv"
abund_dat <- read.csv(infile, stringsAsFactors = F)

s_pie_large <- frag_div %>% 
   filter(S_PIE_mean > S_obs)
plot(S_PIE_mean ~ S_obs, data = frag_div, log = "xy")
abline(0,1)

abund_s_pie_large <- abund_dat %>%
   inner_join(s_pie_large) %>%
   select(dataset_label, frag_id, sample_id, species, abundance, N, S_obs, S_PIE_mean)
dim(abund_s_pie_large)
View(abund_s_pie_large)



# Get numbers of fragments per study
 frag_div %>%
   select(dataset_label, frag_id) %>%
   distinct() %>% 
   group_by(dataset_label) %>%
   count() %>% 
   arrange(desc(n))


summary(frag_div)

plot(S_std1_mean ~ S_std2_mean, data = frag_div)
cor(frag_div$S_std1_mean, frag_div$S_std2_mean)
cor_mat <- cor(select(frag_div, S_obs, S_std1_mean:S_chao_mean), use = "pairwise")
write.csv(cor_mat,"correlation_matrix.csv")
pairs(select(frag_div, S_obs, S_std1_mean:S_chao_mean))

# get number of studies with specific sampling design
frag_div %>% 
   select(dataset_label, sample_design) %>%
   distinct() %>% count(sample_design)

sample_designs <- frag_div %>% 
   select(dataset_label, sample_design) %>%
   distinct() %>% 
   arrange(sample_design, dataset_label)

outfile <- path2Dropbox %+% "files_datapaper/Long_format_database/sample_designs_overview.csv"
write.csv(sample_designs, outfile, row.names = F)

# get numbers of samples per fragmennt

sample_designs <- frag_div %>% 
   select(dataset_label, sample_design) %>%
   distinct() %>% 
   arrange(sample_design, dataset_label)
   


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


# compare results with different methods
infile1 <- path2Dropbox %+% "/files_datapaper/Analysis/1_biodiv_frag_fcont_10_mabund_as_is.csv"
frag_div1 <- read.csv(infile1, stringsAsFactors = F)

infile2 <- path2Dropbox %+% "/files_datapaper/Analysis_June_26/2_biodiv_frag_fcont_10_mabund_as_is.csv"
frag_div2 <- read.csv(infile2, stringsAsFactors = F)

frag_div1 <- frag_div1 %>%
   arrange(dataset_label, frag_id)

frag_div2 <- frag_div2 %>%
   arrange(dataset_label, frag_id)


names(frag_div1)
names(frag_div2)

plot(frag_div1$N_std ~ frag_div2$N_std)
plot(frag_div1$S_std ~ frag_div2$S_std_1)
plot(frag_div1$S_std ~ frag_div2$S_std_2)
plot(frag_div1$S_obs ~ frag_div2$S_obs)
abline(0,1)

plot(frag_div1$S_n ~ frag_div2$S_n)
abline(0,1)

plot(frag_div1$S_PIE ~ frag_div2$S_PIE)
abline(0,1)
