############################################################################
### set directories
############################################################################

.setwdntemp <- function(){
   cu <- Sys.info()["user"]
   cn <- Sys.info()["nodename"]
   
   if (cu == "kg83hyby") #KG 
   {
      path2wd <- "C:/Users/kg83hyby/Documents/GitHub/FragFrame_1/" 
      path2Dropbox <- "C:/Users/kg83hyby/Dropbox/Habitat loss meta-analysis/"
      path2temp <- "C:/Users/kg83hyby/Documents/temp/FragFrame_1/" 
   }  
   else {#FM
      path2wd <- "C:/Users/May/Documents/FelixMay/Fragmentation_Extinction/FragFrame_1/" 
      path2Dropbox <- "C:/Users/May/Dropbox (iDiv)/Habitat loss meta-analysis/"
      path2temp <- "C:/Users/May/Dropbox (iDiv)/Habitat loss meta-analysis/analysis/"
   }
   return(list(path2temp,path2Dropbox,path2wd))
}

set.list <-  .setwdntemp()
path2temp <- set.list[[1]]
path2Dropbox <- set.list[[2]]
path2wd <- set.list[[3]]

############################################################################
### some helper functions
############################################################################
### helper function to combine strings
"%+%" <- function(x,y) paste(x,y,sep="")

### helper function to try things out
is.error <- function(x) inherits(x, "try-error")

############################################################################
### Load and install all needed libraries
############################################################################
needed_libs <- c("devtools", # download from github
                 "mobr", # calculation of biodiversity indices
                 "vegan", # for diversity indices
                 "adespatial", # for beta-diversity partitioning
                 "car", # for logit transformation
                 "metafor", # for meta-analysis 
                 "lme4", # for lmer
                 "ggplot2", # for plotting
                 "gridExtra", # for multiple plots using grid.arrange()
                 "grid", # for extracting legends with grid.draw()
                 "plot3D", # for multidimensional plotting
                 "plyr",
                 "dplyr", # for data manipulation
                 "reshape2", # for restructuring datasets with melt()
                 "xlsx",  # for reading Excel spreadsheets
                 "SpadeR",
                 "RColorBrewer",
                 "raster"# for trim()
)

usePackage <- function(p) {
   #if(p=="iNEXT")   install_github('JohnsonHsieh/iNEXT') # iNEXT is on CRAN !
   if (!is.element(p, installed.packages()[,1])) {   
      if(p == "mobr") {install_github('MoBiodiv/mobr')}  
      if(p == "adespatial") {install.packages("XML",type="binary")}
      install.packages(p, dep = TRUE)
   }
   require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)

### document system properties
session_info()

rm(needed_libs, usePackage)

