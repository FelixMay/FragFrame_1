div_df <- read.csv(path2temp %+% "DiversityData.csv", sep=",")
#str(div_df)

# div_df <- subset(div_df, base_cov>0.5) # exclude studies with base coverage less than 0.5 that indicates incomplete sampling
# We do not calculate a base coverage anymore, because we (Jon + Felix) decided to use
# Schao instead of coverage standardized D0hat (for consistency with mob framework)
# So we might need another quality check

# BDmetrics <- c("S_obs","D0_hat","N_std","ENS_pie")
BDmetrics <- c("N", "N_std","S_obs","S_std","S_n1","S_n2","S_asymp","S_PIE")


###############################################################
### 1. largest vs smallest fragment incl continuous as fragment using log RR
###############################################################
#ES_frag_df <- data.frame(Case.ID = unique(div_df$filename))
ES_frag_df <- unique(dplyr::select(div_df, filename, sample_design,
                                   repl_part_BS_qF, repl_part_BS_qT)
                     )

names(ES_frag_df)[1] <- "Case.ID"
namES_df <- c("n.fragment","log.rr.entity.size",c("ES.") %+% BDmetrics)
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
#ES_frag_group_df <- data.frame(Case.ID=unique(div_df$filename))
ES_frag_group_df <- unique(dplyr::select(div_df, filename, sample_design,
                                         repl_part_BS_qF, repl_part_BS_qT)
                          )

names(ES_frag_group_df)[1] <- "Case.ID"
namES_df <- c("n.fragment","log.rr.entity.size",c("ES.") %+% BDmetrics)
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
#ES_df <- data.frame(Case.ID=unique(div_df$filename))
ES_df <- unique(dplyr::select(div_df, filename, sample_design,
                              repl_part_BS_qF, repl_part_BS_qT)
               )

names(ES_df)[1] <- "Case.ID"
namES_df <- c("n.fragment",c("ES.","ES.var.") %+% rep(BDmetrics,each=2))
ES_df[,namES_df] <- NA

for(i in 1:length(ES_df$Case.ID)){
   print(ES_df$Case.ID[i], max.levels=0)
   sub.df <- subset(div_df, filename==ES_df$Case.ID[i])
   ES_df[i,"n.fragment"] <- nrow(sub.df)

   for(j in BDmetrics){
      r <- ifelse(ES_df[i,"n.fragment"] < 90,
                  2*sin(pi*cor(sub.df$entity.size.rank,sub.df[,j],method="spearman", use = "na.or.complete")/6),
                  cor(sub.df$entity.size.rank,sub.df[,j],method="spearman", use = "na.or.complete")) ### transform rank-correlation into Pearson-correlation, cf. Box 13.3. p 201 in Koricheva et al 2013
      ES_df[i,"ES." %+% j] <- 1/2*log((1+r)/(1-r))
      ES_df[i,"ES.var." %+% j] <- ifelse(ES_df[i,"n.fragment"]>3,1/(ES_df[i,"n.fragment"]-3),NA)
   }
   
   BDmetrics <- c("N", "N_std","S_obs","S_std","S_n1","S_n2","S_asymp","S_PIE")
   fig_name <- path2temp %+% "BoxPlots/" %+% paste(unique(sub.df$filename), ".png", sep = "")
   png(fig_name, width = 10, height = 6, units = "in", res = 200)
   op <- par(mfrow = c(2,5), oma = c(0,0,4,0), las = 1)
   boxplot(N ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "N")
   boxplot(sample_effort ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "Sampling effort")
   boxplot(N_std ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "N standardized")
   boxplot(S_obs ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "Observed S")
   if(!all(is.na(sub.df$S_std))){
      boxplot(S_std ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "Rarefied S")
   }
   if(!all(is.na(sub.df$S_n1))){
      boxplot(S_n1 ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
              main = "Rarefied S to N_min")
   }
   if(!all(is.na(sub.df$S_n2))){
      boxplot(S_n2 ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "Rarefied S to 2*N_min")
   }
   boxplot(S_asymp ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "Asymptotic S")
   boxplot(S_PIE ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "S_PIE")
   if (ES_df$sample_design[i] != "pooled")
      try(boxplot(repl_part_BS_qT ~ entity.size.rank, data = sub.df, xlab = "Fragment size rank",
           main = "repl_part_BS_qT"))
   
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

# Check correlation among BD indices
pairs(ES_frag_group_df[, c("ES.N", "ES.N_std","ES.S_obs","ES.S_std","ES.S_n1","ES.S_n2","ES.S_asymp","ES.S_PIE")])
cor(ES_frag_group_df[, c("ES.N", "ES.N_std","ES.S_obs","ES.S_std","ES.S_n1","ES.S_n2","ES.S_asymp","ES.S_PIE")],
    use = "pairwise.complete.obs")

#pairs(ES_frag_group_df[, c("ES.N_std","ES.S_obs","ES.D0_hat","ES.ENS_pie")])
# cor(ES_frag_group_df[, c("ES.N_std","ES.D0_hat","ES.S_obs","ES.ENS_pie")], use = "pairwise.complete.obs")

# Check correlation among turnover indices
pairs(ES_frag_group_df[, c("repl_part_BS_qT","repl_part_BS_qF")])
cor(ES_frag_group_df[, c("repl_part_BS_qT","repl_part_BS_qF")], use = "pairwise.complete.obs")
