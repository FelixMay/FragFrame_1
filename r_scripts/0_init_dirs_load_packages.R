############################################################################
### set directories and load packages
############################################################################

# Set user dependent working directories

# This is the ONE AND ONLY PLACE where an absolute path should be used!!!

# Shane and Alban, add your working directories here accordingly

user <- Sys.info()["user"]

path2wd <- switch(user,
                  "May" = "C:/Users/May/Documents/FelixMay/Fragmentation_Extinction/FragFrame_1",
                  "sb25gaqy" = "~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/"
)

setwd(path2wd)
   
############################################################################
### some helper functions
############################################################################
### helper function to combine strings

# Are these needed??? I do not think so!!!
#"%+%" <- function(x,y) paste(x, y, sep="")

### helper function to try things out
#is.error <- function(x) inherits(x, "try-error")

############################################################################
### Load all needed libraries
############################################################################
needed_libs <- c("adespatial" , # for beta-diversity partitioning
                 "brms",
                 "cowplot",
                 "devtools",
                 "lme4", # for lmer 
                 "mobr", # calculation of biodiversity indices,
                 "mobsim",
                 "tidyverse",
                 "stringr"
                 #"vegan", # for diversity indices
                 #"iNEXT", # for coverage standardized richness
)

usePackage <- function(p) {
   if (!is.element(p, installed.packages()[,1])) {   
      if(p == "mobr") {install_github('MoBiodiv/mobr')}  
      install.packages(p, dep = TRUE)
   }
   require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)

### document system properties
session_info() 

