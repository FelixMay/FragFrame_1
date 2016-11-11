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
      path2temp <- "C:/Users/kg83hyby/Documents/temp/FragFrame_1" 
   }  
   else {#FM
      path2wd <-"c:/dropbox/fm28towy/Dropbox/Habitat loss meta-analysis/good_datasets/" 
      path2Dropbox <- "c:/dropbox/fm28towy/Dropbox/Habitat loss meta-analysis/"
      path2temp <- ""
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
"%+%" <- function(x,y)paste(x,y,sep="")

### helper function to try things out
is.error <- function(x) inherits(x, "try-error")

############################################################################
### Load and install all needed libraries
############################################################################
needed_libs <- c("devtools", # download from github
                 "MoBspatial", # simulation of species communities
                 "iNEXT", # computes diversity estimates for rarefied and extrapolated samples
                 "vegan", # for diversity indices
                 "ggplot2", # for plotting
                 "gridExtra", # for multiple plots using grid.arrange()
                 "dplyr" # for data manipulation
)
usePackage <- function(p) {
   if(p == "MoBspatial")    install_github('MoBiodiv/MoBspatial')    
   if(p=="iNEXT")   install_github('JohnsonHsieh/iNEXT')
   if (!is.element(p, installed.packages()[,1]))    install.packages(p, dep = TRUE)
   require(p, character.only = TRUE)
}
sapply(needed_libs,usePackage)

rm(needed_libs, usePackage)
############################################################################
### DATA PREPARATION
############################################################################
source(path2wd %+% "Summary_all_files.r")

############################################################################
###  DATA ANALYSIS
############################################################################


############################################################################
###  Plotting
############################################################################
