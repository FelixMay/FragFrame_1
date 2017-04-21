div_df <- read.csv(path2temp %+% "DiversityData.csv", sep=",")
ES_frag_df <- read.csv(file=path2temp %+% "ES_frag_df.csv")
ES_frag_group_df <- read.csv(file=path2temp %+% "ES_frag_group_df.csv")
ES_df <- read.csv(file=path2temp %+% "ES_df.csv")
meta_df <- read.csv(file=path2temp %+% "metaData.csv")

BDmetrics <- c("N_std","S","D0_hat","ENS_pie")

############################################################################
### 1. subset original dataset
############################################################################
ES_frag_df <- ES_frag_df[,c("Case.ID","Study.ID", "n.fragment",
                            sapply(BDmetrics, function(x) paste(c("ES."),x,sep="")))]

ES_frag_group_df <- ES_frag_group_df[,c("Case.ID","Study.ID", "n.fragment",
                            sapply(BDmetrics, function(x) paste(c("ES."),x,sep="")))]

ES_df <- subset(ES_df, n.fragment>3)
ES_df <- ES_df[,c("Case.ID","Study.ID", "n.fragment",
                  sapply(BDmetrics, function(x) paste(c("ES.","ES.var."),x,sep="")))]

############################################################################
### 2.append meta-data
############################################################################
ES_frag_df.complete <- left_join(ES_frag_df,meta_df[,-1],by="Case.ID")
ES_frag_group_df.complete <- left_join(ES_frag_group_df,meta_df[,-1],by="Case.ID")
ES_df.complete <- left_join(ES_df,meta_df[,-1],by="Case.ID")

# ############################################################################
# ### 3. drop unused levels
# ############################################################################
# ES_df.complete$taxa <- factor(ES_df.complete$taxa)[drop=T]
# ES_df.complete$country <- factor(ES_df.complete$country)[drop=T]
# ES_df.complete$continent <- factor(ES_df.complete$continent)[drop=T]
# ES_df.complete$biome <- factor(ES_df.complete$biome)[drop=T]

############################################################################
### 4. restructure datasets _df.complete
############################################################################
## transform into long dataframes
ES_frag_df.complete_long <- melt(ES_frag_df.complete,variable.name="ES",measure.vars=sapply(BDmetrics, function(x) paste(c("ES."),x,sep="")))

ES_frag_group_df.complete_long <- melt(ES_frag_group_df.complete,variable.name="ES",measure.vars=sapply(BDmetrics, function(x) paste(c("ES."),x,sep="")))

ES_df.complete_long <- melt(ES_df.complete,variable.name="ES",measure.vars=sapply(BDmetrics, function(x) paste("ES.",x,sep="")))
ES_df.complete_long$ES.var <- NA
for (BD in BDmetrics){
   ES_df.complete_long$ES.var[ES_df.complete_long$ES=="ES." %+% BD] <- ES_df.complete_long[ES_df.complete_long$ES=="ES." %+% BD,"ES.var." %+% BD]
}
ES_df.complete_long <- ES_df.complete_long[,-which(names(ES_df.complete_long) %in% c("ES.var." %+% BDmetrics))]

############################################################################
### 5. save output
############################################################################

write.csv(ES_frag_df.complete, file=path2temp %+% "ES_frag_df.complete.csv")
write.csv(ES_frag_group_df.complete, file=path2temp %+% "ES_frag_group_df.complete.csv")
write.csv(ES_df.complete, file=path2temp %+% "ES_df.complete.csv")



