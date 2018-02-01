############################################################################
### 00 Initialize Directories & Load Packages And Functions
############################################################################
source("00_InitializeDirectories_LoadPackagesAndFunctions.R", echo=TRUE)

sink(file=path2temp %+% "dataprep_documentation.txt")
############################################################################
### 01 DATA PREPARATION
############################################################################
# 1. create table with raw data Case.ID-Fragment-Rank-Abundance
source(path2wd %+% "01.1_CalculateBDFromAbundance.R") 

# 2. create table with Case.ID - meta-data
source(path2wd %+% "01.2_AddMetaData.r") 

# 3. calculate effect sizes per Case.ID
source(path2wd %+% "01.3_CalculateEffectSizes.r") 

# 4. merge dataframes with meta-data
rm(list=ls())
source("00_InitializeDirectories_LoadPackagesAndFunctions.R", echo=F)
source(path2wd %+% "02.1_DataPrep4Analysis.r") 
save.image(file=path2temp %+% "02.1_Data4Analysis_out.Rdata")

# 5. 
rm(list=ls())
source("00_InitializeDirectories_LoadPackagesAndFunctions.R", echo=F)
source(path2wd %+% "02.2_DescriptiveStats.r") 

sink()

############################################################################
###  02 DATA ANALYSIS
############################################################################
# 1. 
rm(list=ls())
source("00_InitializeDirectories_LoadPackagesAndFunctions.R", echo=F)
#source(path2wd %+% "02.3a_BDAnalysis.r") 
#source(path2wd %+% "02.3b_BetaDivAnalysis.r")
source(path2wd %+% "02.3_DataAnalysis.r") 

save.image(file=path2temp %+% "02.3_DataAnalysis_out.Rdata")

# 2. 
rm(list=ls())
source("00_InitializeDirectories_LoadPackagesAndFunctions.R", echo=F)
sink(file=path2temp %+% "CheckResults/ModelDiagnostics_Output.txt")
source(path2wd %+% "02.4_CheckResults.r") 
sink()

rm(list=ls())
source("00_InitializeDirectories_LoadPackagesAndFunctions.R", echo=F)
source(path2wd %+% "02.5_VisualizeResults.r") 

rm(list=ls())
source("00_InitializeDirectories_LoadPackagesAndFunctions.R", echo=F)
source(path2wd %+% "02.6_SensitivityAnalysis.r") 

