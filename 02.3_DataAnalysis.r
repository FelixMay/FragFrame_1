load(path2temp %+% "Data4Analysis_out.Rdata") 
ls()

############################################################################
### 1. set reference levels to the most common levels
############################################################################
setRefToMostCommonLevel <- function(f) {
   f <- as.factor(f)
   t <- table(f)
   f.new <- relevel(f,ref=as.integer(which(t>=max(t))[[1]]))
   return(f.new)
}

############################################################################
### 2. Analysis pipeline
############################################################################
analysis_func <- function(df, covar="intercept.only", method){#,df_long
 
   ### set reference levels
   for(col in c("taxa","country", "continent", "biome", "fragment.biome","matrix.biome", "fragment.veg", "matrix.veg")){
      df[,col] <- setRefToMostCommonLevel(df[,col])
   }
   
   model <- list()

   if(covar=="intercept.only"){
      mods.formula <- "~1"
   }
   if(covar %in% c("taxa", "country", "continent", "biome", "fragment.biome","matrix.biome", "fragment.veg", "matrix.veg" )){
      mods.formula <- "~" %+% covar %+% "-1"
   }
   
   ############################################################################
   ### for unweighted meta-analysis
   ### Weighted estimation (with inverse-variance weights) is used by default. User-defined weights can be supplied via the weights argument. Unweighted estimation can be used by setting weighted=FALSE in rma.uni. This is the same as setting the weights equal to a constant (Source: R Help rma.uni)
   if(!any(grepl("var",names(df)))){
      df <- data.frame(df,
                       ES.var.S = rep(1,nrow(df)),
                       ES.var.D0_hat = rep(1,nrow(df)),
                       ES.var.N_std = rep(1,nrow(df)),
                       ES.var.ENS_pie = rep(1,nrow(df)))
   }      
      # 1.	Species loss: Fragmentation leads to decrease of observed species richness, i.e. naïve species richness without considering methodological differences.  
      model[["Naive species loss"]] <- rma.mv(yi=ES.S, V=ES.var.S, mods = as.formula(mods.formula), random = ~ 1 | Study.ID, struct="UN", data=df, method=method)
      
      # 2.	Fragmentation per se:
      #    a.	Species loss due to less coverage standardized species richness: Fragmentation leads to decrease in coverage standardized species richness.
      model[["Coverage standardized species loss"]] <- rma.mv(yi=ES.D0_hat, V=ES.var.D0_hat, mods = as.formula(mods.formula), random = ~ 1 | Study.ID, struct="UN", data=df, method=method)
      
      #     b.	Species loss due to less Individuals: Fragmentation decreases the number of individuals.
      model[["Less individuals"]] <- rma.mv(yi=ES.N_std,V=ES.var.N_std, mods = as.formula(mods.formula), random = ~ 1 | Study.ID, struct="UN", data=df, method=method)
      #     c.	Species loss due to PIE/ENS_PIE: Fragmentation reduces the evenness.
      model[["Lower evenness"]] <- rma.mv(yi=ES.ENS_pie,V=ES.var.ENS_pie, mods = as.formula(mods.formula), random = ~ 1 | Study.ID, struct="UN", data=df, method=method)
      
      # # 3.	Common vs rare: Fragmentation leads to decrease in coverage standardized species richness and to a change in the species-abundance distribution, e.g. rare species are more affected than common species (increasing ENS_PIE).
      # df_long_sub <- subset(df_long, ES %in% c("ES.D0_hat", "ES.ENS_pie"))
      # df_long_sub$ES <- factor(df_long_sub$ES)
      # model[["Common vs. rare"]] <- rma.mv(yi=value,V=ES.var, mods = as.formula(mods.formula.extended), random = ~ ES | Study.ID, struct="UN", data=df_long_sub, method=method)
   return(model)
}

# 1. Grand Mean
model_frag <- analysis_func(df=ES_frag_df.complete, method="REML")
model_frag_group <- analysis_func(df=ES_frag_group_df.complete, method="REML")
model_gradient <- analysis_func(df=ES_df.complete, method="REML")

# 2.	Taxonomic group
#      Do effects vary systematically vary with taxonomic group, i.e. is the overall heterogeneity reduced within groups of the same taxa?
modelTaxa_frag <- analysis_func(df=ES_frag_df.complete, covar="taxa", method="REML")
modelTaxa_frag_group <- analysis_func(df=ES_frag_group_df.complete,covar="taxa", method="REML")
modelTaxa_gradient <- analysis_func(df=ES_df.complete, covar="taxa", method="REML")

# 3.	Habitat and biome of fragments and matrix
#     Do effects vary systematically vary with types of habitat and biome of fragments and surrounding matrix, i.e. is the overall heterogeneity reduced within groups of the same habitat/biome type of fragment/matrix or combinations?
### TO DO: Come up with useful categories of habitat and biome of fragments and matrix


