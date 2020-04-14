############################################################################
### set directories and load packages
############################################################################

.setwdntemp <- function(){
   cu <- Sys.info()["user"]
   cn <- Sys.info()["nodename"]
   
   # if (cu == "kg83hyby") #KG 
   # {
   #    path2wd <- "C:/Users/kg83hyby/Documents/GitHub/FragFrame_1/" 
   #    path2Dropbox <- "C:/Users/kg83hyby/Dropbox/Habitat loss meta-analysis/"
   #    path2temp <- "C:/Users/kg83hyby/Documents/temp/FragFrame_1/" 
   # }  
   if (cu == 'sb25gaqy')
   {
      path2wd <- '~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/'
      path2Dropbox <- '~/Dropbox/Frag Database (new)'
      path2data <- '~/Dropbox/Frag Database (new)/files_datapaper/Analysis/'
      path2meta <- '~/Dropbox/Frag Database (new)/'
      path2temp <- '/Users/sb25gaqy/Dropbox/1current/fragmentation_synthesis/temp/'
   }
   else {#FM
      path2wd <- "C:/Users/May/Documents/FelixMay/Fragmentation_Extinction/FragFrame_1/" 
      path2Dropbox <- "C:/Users/May/Dropbox (Privat)/Frag Database (new)/"
      path2data <- "C:/Users/May/Dropbox (Privat)/Frag Database (new)/files_datapaper/Analysis/"
      path2meta <- "C:/Users/May/Dropbox (Privat)/Frag Database (new)/"
      path2temp <- "C:/Users/May/Dropbox (Privat)/Frag Database (new)/Analysis/"
   }
   return(list(path2temp,path2Dropbox,path2wd, path2data, path2meta))
}

set.list <-  .setwdntemp()
path2temp <- set.list[[1]]
path2Dropbox <- set.list[[2]]
path2wd <- set.list[[3]]
path2data <- set.list[[4]]
path2meta <- set.list[[5]]

############################################################################
### some helper functions
############################################################################
### helper function to combine strings
"%+%" <- function(x,y) paste(x, y, sep="")

### helper function to try things out
is.error <- function(x) inherits(x, "try-error")

############################################################################
### Load and install all needed libraries
############################################################################
needed_libs <- c(#"devtools", # download from github
                 "mobr", # calculation of biodiversity indices
                 "vegan", # for diversity indices
                 "adespatial" , # for beta-diversity partitioning
                 #"car", # for logit transformation
                 #"metafor", # for meta-analysis 
                 "lme4", # for lmer
                 #"gridExtra", # for multiple plots using grid.arrange()
                 #"grid", # for extracting legends with grid.draw()
                 #"plot3D", # for multidimensional plotting
                 "tidyverse",
                 #"reshape2", # for restructuring datasets with melt()
                 #"xlsx",  # for reading Excel spreadsheets
                 #"SpadeR",
                 #"RColorBrewer",
                 #"raster", # for trim()
                 "iNEXT", # for coverage standardized richness
                 "brms",
                 "cowplot",
                 "stringr"

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
# session_info() # devtools function maybe not needed here

rm(needed_libs, usePackage)

