div_df <- read.csv(path2temp %+% "DiversityData.csv", sep=";")
ES_df <- read.csv(file=path2temp %+% "ES_df.csv")
meta_df <- read.csv(file=path2temp %+% "metaData.csv")

ES_df <- subset(ES_df, n.fragment>3)
ES_df <- ES_df[,c("Study.ID", "n.fragment",
                  "z.N", "z.var.N", 
                  "z.S", "z.var.S",
                  "z.Shannon", "z.var.Shannon",
                  "z.PIE", "z.var.PIE",
                  "z.S_cov", "z.var.S_cov",
                  "z.Shannon_cov", "z.var.Shannon_cov",
                  "z.PIE_cov", "z.var.PIE_cov")]

## restructure ES_df
df.list <- list()
for(BD in c("N","S","Shannon", "PIE", "S_cov","Shannon_cov", "PIE_cov")){
   df.list[[BD]] <- ES_df[,c("Study.ID", "n.fragment","z." %+% BD, "z.var." %+% BD)]
   names(df.list[[BD]]) <- c("Study.ID", "n.fragment","z", "z.var")
   df.list[[BD]]$BD <- BD
}

df <- bind_rows(df.list)
df <- df[order(df$Study.ID),]
df$BD <- as.factor(df$BD)

## append meta data
df.complete <- left_join(df,meta_df[,-1],by="Study.ID")

## set reference levels to the most common levels
setRefToMostCommonLevel <- function(f) {
   f <- as.factor(f)
   t <- table(f)
   relevel(f,ref=as.integer(which(t>=max(t))[[1]]))
}
df.complete[,c("taxa","country", "continent", "biome")] <- sapply(df.complete[,c("taxa","country", "continent", "biome")],setRefToMostCommonLevel)

## multivariate model, cf. http://www.metafor-project.org/doku.php/analyses:vanhouwelingen2002#bivariate_approach
model <- list()
model[["BD"]] <- rma.mv(z, z.var, mods = ~ BD - 1, random = ~ BD | Study.ID, struct="UN", data=df.complete, method="ML")
model[["Full"]] <- rma.mv(z, z.var, mods = ~ BD + taxa + continent, random = ~ BD | Study.ID, struct="UN", data=df.complete, method="ML")

write.csv(df.complete, file=path2temp %+% "df.complete.csv")
save(df.complete, model, file="MA.model.Rdata")
