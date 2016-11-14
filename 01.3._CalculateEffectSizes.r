div.df <- read.csv(path2temp %+% "DiversityData.csv", sep=";")
# div.df <- read.csv(path2Dropbox %+% "Analysis/DiversityData.csv",sep=";")

#table(div.df$filename)
# table(div.df$entity.id)

### get rid of irrelevant data, e.g. matrix, clearcut
div.df <- subset(div.df, entity.size.rank >0)

### calculate rank-correlation
ES.df <- data.frame(Study.ID=unique(div.df$filename))
names.df <- c("z.","z.var.") %+% names(div.df)[-(1:3)]
ES.df[,names.df] <- NA

for(i in 1:length(ES.df$Study.ID)){
   sub.df <- subset(div.df, filename==ES.df$Study.ID[i])
   for(j in names(div.df)[-(1:3)]){
      r <- ifelse(nrow(sub.df)<90,2*sin(pi*cor(sub.df$entity.size.rank,sub.df[,j])/6),cor(sub.df$entity.size.rank,sub.df[,j])) ### transform rank-correlation into Pearson-correlation, cf. Box 13.3. p 201 in Koricheva et al 2013
      ES.df[i,"z." %+% j] <- 1/2*log((1+r)/(1-r))
      ES.df[i,"z.var." %+% j] <- ifelse(nrow(sub.df)>3,1/(nrow(sub.df)-3),NA)
   }
}

write.csv(ES.df, file=path2temp %+% "ES.df.csv")


### forest plots
forest(x=ES.df$z.N,vi=ES.df$z.var.N,xlim=c(-10,10))
