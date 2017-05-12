div_df <- read.csv(path2temp %+% "DiversityData.csv", sep=",")

str(div_df)

BDmetrics <- c("N_std","S","D0_hat","ENS_pie")

###############################################################
### 1. largest vs smallest fragment incl continuous as fragment using log RR
###############################################################
ES_frag_df <- data.frame(Case.ID=unique(div_df$filename))
namES_df <- c("sample.design","n.fragment","log.rr.entity.size",c("ES.") %+% BDmetrics)
ES_frag_df[,namES_df] <- NA

for(i in 1:length(ES_frag_df$Case.ID)){
   sub.df <- subset(div_df, filename %in% ES_frag_df$Case.ID[i])
   small.df <- sub.df[which(sub.df$entity.size.rank == min(sub.df$entity.size.rank)),]
   large.df <- sub.df[which(sub.df$entity.size.rank == max(sub.df$entity.size.rank)),]
   ES_frag_df[i,"n.fragment"] <- nrow(sub.df)
   ES_frag_df[i,"log.rr.entity.size"] <-
      log(mean(small.df$entity.size)/mean(large.df$entity.size)) #log-ratio of smallest vs. largest fragment area
   n.small <- nrow(small.df)
   n.large <- nrow(large.df) 
   for(j in BDmetrics){
      ES_frag_df[i,"ES." %+% j] <- log(mean(small.df[,j])/mean(large.df[,j]))
   }
}

###############################################################
### 2. largest (incl. continuous) vs smallest fragment group using log RR
###############################################################
ES_frag_group_df <- data.frame(Case.ID=unique(div_df$filename))
namES_df <- c("sample.design","n.fragment","log.rr.entity.size",c("ES.") %+% BDmetrics)
ES_frag_group_df[,namES_df] <- NA

for(i in 1:length(ES_frag_group_df$Case.ID)){
   sub.df <- subset(div_df, filename %in% ES_frag_group_df$Case.ID[i])
   small.df <- sub.df[which(sub.df$entity.size.rank < quantile(unique(sub.df$entity.size.rank),probs=0.25)),]
   large.df <- sub.df[which(sub.df$entity.size.rank > quantile(unique(sub.df$entity.size.rank),probs=0.75)),]
   ES_frag_group_df[i,"n.fragment"] <- nrow(sub.df)
   ES_frag_group_df[i,"log.rr.entity.size"] <-
      log(mean(small.df$entity.size)/mean(large.df$entity.size)) #log-ratio of smallest vs. largest fragment area
   n.small <- nrow(small.df)
   n.large <- nrow(large.df) 
   
   for(j in BDmetrics){
      ES_frag_group_df[i,"ES." %+% j] <- log(mean(small.df[,j])/mean(large.df[,j]))
   }
}

###############################################################
### 3. Gradient of habitat fragmentation using Fishers' z
###############################################################
#calculate rank-correlation
ES_df <- data.frame(Case.ID=unique(div_df$filename))
namES_df <- c("sample.design","n.fragment",c("ES.","ES.var.") %+% rep(BDmetrics,each=2))
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
   
   fig_name <- path2temp %+% "BoxPlots/" %+% paste(unique(sub.df$filename), ".png", sep = "")
   png(fig_name, width = 8, height = 6, units = "in", res = 200)
   op <- par(mfrow = c(2,3), oma = c(0,0,4,0), las = 1)
   boxplot(N ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "N")
   boxplot(sample_effort ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "Sampling effort")
   boxplot(N_std ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "N standardized")
   boxplot(S ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "S")
   boxplot(D0_hat ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "Asymptotic S")
   boxplot(ENS_pie ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "ENS_PIE")
   mtext(unique(sub.df$filename), side = 3,
         line = 2, outer = T, cex = 1.2)
   mtext(paste(ES_df$n.fragment[i]," Fragments, ", ES_df$sample.design[i], sep = ""),
               side = 3, line = 0, outer = T)
   
   dev.off()
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


