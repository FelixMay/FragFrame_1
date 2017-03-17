div_df <- read.csv(path2temp %+% "DiversityData.csv", sep=",")

str(div_df)

BDmetrics <- c("N_std","S","D0_hat","ENS_pie")

###############################################################
### 1. largest vs smallest fragment incl continuous as fragment using log RR
###############################################################
ES_frag_df <- data.frame(Case.ID=unique(div_df$filename))
namES_df <- c("n.fragment",c("ES.") %+% BDmetrics)
ES_frag_df[,namES_df] <- NA

for(i in 1:length(ES_frag_df$Case.ID)){
   sub.df <- subset(div_df, filename %in% ES_frag_df$Case.ID[i])
   small.df <- sub.df[which.min(sub.df$entity.size.rank),]
   large.df <- sub.df[which.max(sub.df$entity.size.rank),]
   ES_frag_df[i,"n.fragment"] <- nrow(sub.df)
   for(j in BDmetrics){
      ES_frag_df[i,"ES." %+% j] <- log(small.df[,j]/large.df[,j])
   }
}

###############################################################
### 2. largest (incl. continuous) vs smallest fragment group using log RR
###############################################################
ES_frag_group_df <- data.frame(Case.ID=unique(div_df$filename))
namES_df <- c("n.fragment",c("ES.") %+% BDmetrics)
ES_frag_group_df[,namES_df] <- NA

for(i in 1:length(ES_frag_group_df$Case.ID)){
   sub.df <- subset(div_df, filename %in% ES_frag_group_df$Case.ID[i])
   small.df <- sub.df[which(sub.df$entity.size.rank<median(sub.df$entity.size.rank)),]
   large.df <- sub.df[which(sub.df$entity.size.rank>median(sub.df$entity.size.rank)),]
   ES_frag_group_df[i,"n.fragment"] <- nrow(sub.df)
   for(j in BDmetrics){
      ES_frag_group_df[i,"ES." %+% j] <- log(mean(small.df[,j])/mean(large.df[,j]))
   }
}

###############################################################
### 3. Gradient of habitat fragmentation using Fishers' z
###############################################################
#calculate rank-correlation
ES_df <- data.frame(Case.ID=unique(div_df$filename))
namES_df <- c("n.fragment",c("ES.","ES.var.") %+% rep(BDmetrics,each=2))
ES_df[,namES_df] <- NA

for(i in 1:length(ES_df$Case.ID)){
   sub.df <- subset(div_df, filename==ES_df$Case.ID[i])
   ES_df[i,"n.fragment"] <- nrow(sub.df)
   for(j in BDmetrics){
      r <- ifelse(ES_df[i,"n.fragment"] < 90,
                  2*sin(pi*cor(sub.df$entity.size.rank,sub.df[,j],method="spearman", use = "na.or.complete")/6),
                  cor(sub.df$entity.size.rank,sub.df[,j],method="spearman", use = "na.or.complete")) ### transform rank-correlation into Pearson-correlation, cf. Box 13.3. p 201 in Koricheva et al 2013
      ES_df[i,"ES." %+% j] <- 1/2*log((1+r)/(1-r))
      ES_df[i,"ES.var." %+% j] <- ifelse(ES_df[i,"n.fragment"]>3,1/(ES_df[i,"n.fragment"]-3),NA)
   }
}

###############################################################
### create unique study identifier
###############################################################
splits <- sapply(as.character(ES_df$Case.ID),function(x) strsplit(x,"_")[[1]])
ES_frag_df$Study.ID <- ES_frag_group_df$Study.ID <- ES_df$Study.ID <- sapply(splits, function(x) paste(x[1], x[2], sep="_"))

### write csv
write.csv(ES_frag_df, file=path2temp %+% "ES_frag_df.csv")
write.csv(ES_frag_group_df, file=path2temp %+% "ES_frag_group_df.csv")
write.csv(ES_df, file=path2temp %+% "ES_df.csv")

### Check studies with no effect size
### Not needed for automated execution, just for error checking
# (check <- ES_df[is.na(ES_df$ES.S), ])
# div_df[div_df$filename %in% check$Case.ID,  ]


