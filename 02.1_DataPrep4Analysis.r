div_df <- read.csv(path2temp %+% "DiversityData.csv", sep=";")
ES_df <- read.csv(file=path2temp %+% "ES_df.csv")
meta_df <- read.csv(file=path2temp %+% "metaData.csv")

############################################################################
### 1. subset original dataset
############################################################################
ES_df <- subset(ES_df, n.fragment>3)
ES_df <- ES_df[,c("Study.ID", "n.fragment",
                  "z.N", "z.var.N", 
                  "z.S", "z.var.S",
                  "z.PIE", "z.var.PIE",
                  "z.S_cov", "z.var.S_cov")]

############################################################################
### 2.append meta-data
############################################################################
ES_df.complete <- left_join(ES_df,meta_df[,-1],by="Study.ID")

############################################################################
### 3. drop unused levels
############################################################################
ES_df.complete$taxa <- factor(ES_df.complete$taxa)[drop=T]
ES_df.complete$country <- factor(ES_df.complete$country)[drop=T]
ES_df.complete$continent <- factor(ES_df.complete$continent)[drop=T]
ES_df.complete$biome <- factor(ES_df.complete$biome)[drop=T]

############################################################################
### 4. restructure ES_df.complete
############################################################################
df.list <- list()
for(BD in c("N","S","PIE", "S_cov")){
   df.list[[BD]] <- ES_df.complete[,c("Study.ID", "n.fragment","z." %+% BD, "z.var." %+% BD, "taxa", "country", "continent", "biome")]
   names(df.list[[BD]]) <- c("Study.ID", "n.fragment","z", "z.var", "taxa", "country", "continent", "biome")
   df.list[[BD]]$BD <- BD
}

df.complete <- bind_rows(df.list)
df.complete <- df.complete[order(df.complete$Study.ID),]
df.complete$BD <- as.factor(df.complete$BD)

write.csv(ES_df.complete, file=path2temp %+% "ES_df.complete.csv")
write.csv(df.complete, file=path2temp %+% "df.complete.csv")

save(ES_df.complete,df.complete, file=path2temp %+% "Data4Analysis.Rdata")


