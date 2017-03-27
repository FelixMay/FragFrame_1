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
      path2wd <-"c:/FelixMay/Fragmentation_Extinction/FragFrame_1/" 
      path2Dropbox <- "c:/dropbox/fm28towy/Dropbox/Habitat loss meta-analysis/"
      path2temp <- "c:/dropbox/fm28towy/Dropbox/Habitat loss meta-analysis/Analysis/"
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
                 "MoBspatial", # simulation of species communities
                 "iNEXT", # computes diversity estimates for rarefied and extrapolated samples
                 "vegan", # for diversity indices
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
                 "SpadeR"
)

usePackage <- function(p) {
   #if(p=="iNEXT")   install_github('JohnsonHsieh/iNEXT') # iNEXT is on CRAN !
   if (!is.element(p, installed.packages()[,1])) {   
      if(p == "MoBspatial") {install_github('MoBiodiv/MoBspatial')}  
      install.packages(p, dep = TRUE)
   }
   require(p, character.only = TRUE)
}

### document system properties
version
# _                           
# platform       x86_64-w64-mingw32          
# arch           x86_64                      
# os             mingw32                     
# system         x86_64, mingw32             
# status                                     
# major          3                           
# minor          2.5                         
# year           2016                        
# month          04                          
# day            14                          
# svn rev        70478                       
# language       R                           
# version.string R version 3.2.5 (2016-04-14)
# nickname       Very, Very Secure Dishes    

ip <- as.data.frame(installed.packages()[needed_libs,c(1,3)])
rownames(ip) <- NULL
ip
# Package    Version
# 1    devtools     1.10.0
# 2  MoBspatial 0.0.0.9000
# 3       iNEXT     2.0.12
# 4       vegan      2.4-1
# 5     metafor      1.9-8
# 6     ggplot2      2.1.0
# 7   gridExtra      2.2.1
# 8       dplyr      0.5.0
# 9        xlsx      0.5.7
# 10     SpadeR      0.1.1

sapply(needed_libs, usePackage)

rm(needed_libs, usePackage)

