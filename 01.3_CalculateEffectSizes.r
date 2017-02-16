div_df <- read.csv(path2temp %+% "DiversityData.csv", sep=",")

str(div_df)

### 1. largest vs smallest fragment incl continuous as fragment
ES_frag_df <- data.frame(Case.ID=unique(div_df$filename))
namES_df <- c("n.fragment",c("logRR.") %+% names(div_df)[-(1:3)])
ES_frag_df[,namES_df] <- NA
# small.df <- large.df <- ES_frag_df
# small.df$type <- "small"
# large.df$type <- "large"

for(i in 1:length(ES_frag_df$Case.ID)){
   sub.df <- subset(div_df, filename %in% ES_frag_df$Case.ID[i])
   small.df <- sub.df[which.min(sub.df$entity.size.rank),]
   large.df <- sub.df[which.max(sub.df$entity.size.rank),]
   for(j in names(div_df)[-(1:3)]){
      ES_frag_df[i,"n.fragment"] <- nrow(sub.df)
      ES_frag_df[i,"logRR." %+% j] <- log(small.df[,j]/large.df[,j])
   }
}
#raw_frag_df <- rbind(small.df,large.df) # raw data of small vs large fragments to be analysed in a glmer-framework, cf Newbold et al. 2015

### 2. largest vs smallest fragment group
ES_frag_group_df <- data.frame(Case.ID=unique(div_df$filename))
namES_df <- c("n.fragment",c("logRR.") %+% rep(names(div_df)[-(1:3)],each=2))
ES_frag_group_df[,namES_df] <- NA

for(i in 1:length(ES_frag_group_df$Case.ID)){
   sub.df <- subset(div_df, filename %in% ES_frag_group_df$Case.ID[i])
   small.df <- sub.df[which(sub.df$entity.size.rank<median(sub.df$entity.size.rank)),]
   large.df <- sub.df[which(sub.df$entity.size.rank>median(sub.df$entity.size.rank)),]
   for(j in names(div_df)[-(1:3)]){
      ES_frag_group_df[i,"n.fragment"] <- nrow(sub.df)
      ES_frag_group_df[i,"logRR." %+% j] <- log(mean(small.df[,j])/mean(large.df[,j]))
   }
}

### 3. Gradient of habitat fragmentation
#calculate rank-correlation
ES_df <- data.frame(Case.ID=unique(div_df$filename))
namES_df <- c("n.fragment",c("z.","z.var.") %+% names(div_df)[-(1:3)])
ES_df[,namES_df] <- NA

for(i in 1:length(ES_df$Case.ID)){
   sub.df <- subset(div_df, filename==ES_df$Case.ID[i])
   for(j in names(div_df)[-(1:3)]){
      ES_df[i,"n.fragment"] <- nrow(sub.df)
      r <- ifelse(ES_df[i,"n.fragment"] < 90,
                  2*sin(pi*cor(sub.df$entity.size.rank,sub.df[,j],method="spearman", use = "na.or.complete")/6),
                  cor(sub.df$entity.size.rank,sub.df[,j],method="spearman", use = "na.or.complete")) ### transform rank-correlation into Pearson-correlation, cf. Box 13.3. p 201 in Koricheva et al 2013
      ES_df[i,"z." %+% j] <- 1/2*log((1+r)/(1-r))
      ES_df[i,"z.var." %+% j] <- ifelse(ES_df[i,"n.fragment"]>3,1/(ES_df[i,"n.fragment"]-3),NA)
   }
}

### create unique study identifier
splits <- sapply(as.character(ES_df$Case.ID),function(x) strsplit(x,"_")[[1]])
ES_frag_df$Study.ID <- ES_frag_group_df$Study.ID <- ES_df$Study.ID <- sapply(splits, function(x) paste(x[1], x[2], sep="_"))

### write csv
write.csv(ES_frag_df, file=path2temp %+% "ES_frag_df.csv")
write.csv(ES_frag_group_df, file=path2temp %+% "ES_frag_group_df.csv")
write.csv(ES_df, file=path2temp %+% "ES_df.csv")

### Check studies with no effect size
### Not needed for automated execution, just for error checking
# (check <- ES_df[is.na(ES_df$z.S), ])
# div_df[div_df$filename %in% check$Case.ID,  ]


