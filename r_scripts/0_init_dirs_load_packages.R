############################################################################
### set directories and load packages
############################################################################

# Set user dependent working directories

user <- Sys.info()["user"]
path2wd <- switch(user,
                  "May" = "C:/Users/May/Documents/FelixMay/Fragmentation_Extinction/FragFrame_1/",
                  "sb25gaqy" = "~/Dropbox/1current/fragmentation_synthesis/FragFrame_1/",
						"as80fywe" = "C:/Users/as80fywe/idiv/mob/mobsim/FragFrame_1/"
)

setwd(path2wd)
   
############################################################################
### helper functions
############################################################################

### helper function to try things out
is.error <- function(x) inherits(x, "try-error")

############################################################################
### Load all needed libraries
############################################################################

needed_libs <- c("adespatial" , # for beta-diversity partitioning
                 "brms",
                 "cowplot",
                 "devtools",
                 "iNEXT", # for coverage standardized richness
                 "ggridges", 
                 "mobr", # calculation of biodiversity indices
					       "mobsim", # for random community simulation
                 "tidyverse",
                 "stringr"
                 #"vegan", # for diversity indices
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

