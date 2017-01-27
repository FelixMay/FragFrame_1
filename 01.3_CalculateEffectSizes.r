div_df <- read.csv(path2temp %+% "DiversityData.csv", sep=",")

str(div_df)
### calculate rank-correlation
ES_df <- data.frame(Study.ID=unique(div_df$filename))
namES_df <- c("n.fragment",c("z.","z.var.") %+% rep(names(div_df)[-(1:3)],each=2))
ES_df[,namES_df] <- NA

for(i in 1:length(ES_df$Study.ID)){
   sub.df <- subset(div_df, filename==ES_df$Study.ID[i])
   for(j in names(div_df)[-(1:3)]){
      ES_df[i,"n.fragment"] <- nrow(sub.df)
      r <- ifelse(ES_df[i,"n.fragment"]<90,2*sin(pi*cor(sub.df$entity.size.rank,sub.df[,j])/6),cor(sub.df$entity.size.rank,sub.df[,j])) ### transform rank-correlation into Pearson-correlation, cf. Box 13.3. p 201 in Koricheva et al 2013
      ES_df[i,"z." %+% j] <- 1/2*log((1+r)/(1-r))
      ES_df[i,"z.var." %+% j] <- ifelse(ES_df[i,"n.fragment"]>3,1/(ES_df[i,"n.fragment"]-3),NA)
   }
}

write.csv(ES_df, file=path2temp %+% "ES_df.csv")

