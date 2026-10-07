############################################################################
### set directories and load packages
############################################################################

# Paths are relative to the project root (open FragFrame_1.Rproj)
path2wd <- "./"

############################################################################
### helper functions
############################################################################

### helper function to try things out
is.error <- function(x) inherits(x, "try-error")

############################################################################
### Load all needed libraries
############################################################################

needed_libs <- c(
  "adespatial", # for beta-diversity partitioning
  "brms",
  "cowplot",
  "devtools",
  "iNEXT", # for coverage standardized richness
  "ggridges",
  "ggstance",
  "gridGraphics",
  "magick",
  "maps",
  "mobr", # calculation of biodiversity indices
  "mobsim", # for random community simulation
  "tidyverse",
  "stringr",
  "viridisLite"
)

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[, 1])) {
    if (p == "mobr") {
      install_github('MoBiodiv/mobr')
    }
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)

### document system properties
session_info()
