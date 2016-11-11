div.df <- read.csv(path2temp %+% "DiversityData.csv", sep=";")
# div.df <- read.csv(path2Dropbox %+% "Analysis/DiversityData.csv",sep=";")

#table(div.df$filename)
# table(div.df$entitiy.id)

### get rid of irrelevant data, e.g. matrix, clearcut
div.df <- subset(div.df, entitiy.size.rank >0)

### calculate rank-correlation

### transform rank-correlation into Pearson-correlation, cf. Box 13.3. p 201 in Koricheva et al 2013