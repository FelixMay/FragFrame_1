############################################################################
### 01 DATA PREPARATION
############################################################################
# 1. create table with raw data Study.ID-Fragment-Rank-Abundance
source(path2wd %+% "01.1_CalculateBDFromAbundance.r") 

# 2. create table with Study.ID - meta-data
source(path2wd %+% "01.2_AddMetaData.r") 

# 3. calculate effect sizes per Study.ID and merge with meta-data
source(path2wd %+% "01.3_CalculateEffectSizes.r") 

############################################################################
###  02 DATA ANALYSIS
############################################################################
# 1. Forest Plots of effect sizes and Pairwise Correlation plot of effect sizes
source(path2wd %+% "02.1_DescriptiveStats.r") 

############################################################################
###  03 PLOTTING
############################################################################
