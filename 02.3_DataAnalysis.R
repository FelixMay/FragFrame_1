load(path2temp %+% "02.1_Data4Analysis_out.Rdata") 
ls()

#--------------------------------------------------------
### 1. set reference levels to the most common levels
setRefToMostCommonLevel <- function(f) {
   f <- as.factor(f)
   t <- table(f)
   f.new <- relevel(f,ref=as.integer(which(t>=max(t))[[1]]))
   return(f.new)
}

#--------------------------------------------------------
### 2. Analysis pipeline
analysis_func <- function(df, mods.formula, method){#,df_long
 
   ### set reference levels
   for(col in c("taxa","biome","matrix.category","time.since.fragmentation")){
      df[,col] <- setRefToMostCommonLevel(df[,col])
   }
   
   model <- list()

   ############################################################################
   ### for unweighted meta-analysis
   ### Weighted estimation (with inverse-variance weights) is used by default. User-defined weights can be supplied via the weights argument. Unweighted estimation can be used by setting weighted=FALSE in rma.uni. This is the same as setting the weights equal to a constant (Source: R Help rma.uni)
   if(!any(grepl("var",names(df)))){
      df <- data.frame(df,
                       ES.var.D0_hat = rep(1,nrow(df)),
                       ES.var.N_std = rep(1,nrow(df)),
                       ES.var.ENS_pie = rep(1,nrow(df)))
   }      
      # Fragmentation per se:
      # #    a.	Species loss due to less coverage standardized species richness: Fragmentation leads to decrease in coverage standardized species richness.
      # model[["Coverage standardized species loss"]] <- rma.mv(yi=ES.D0_hat, V=ES.var.D0_hat, mods = as.formula(mods.formula), random = ~ 1 | Study.ID, struct="UN", data=df, method=method)
      # #     b.	Species loss due to less Individuals: Fragmentation decreases the number of individuals.
      # model[["Less individuals"]] <- rma.mv(yi=ES.N_std,V=ES.var.N_std, mods = as.formula(mods.formula), random = ~ 1 | Study.ID, struct="UN", data=df, method=method)
      # #     c.	Species loss due to PIE/ENS_PIE: Fragmentation reduces the evenness.
      # model[["Lower evenness"]] <- rma.mv(yi=ES.ENS_pie,V=ES.var.ENS_pie, mods = as.formula(mods.formula), random = ~ 1 | Study.ID, struct="UN", data=df, method=method)
   
      #    a.	Species loss due to less coverage standardized species richness: Fragmentation leads to decrease in coverage standardized species richness.
      model[["S_asym"]] <- rma.mv(yi=ES.D0_hat, V=ES.var.D0_hat, mods = as.formula(mods.formula), random = ~ 1 | Study.ID, struct="UN", data=df, method=method)
      #     b.	Species loss due to less Individuals: Fragmentation decreases the number of individuals.
      model[["N"]] <- rma.mv(yi=ES.N_std,V=ES.var.N_std, mods = as.formula(mods.formula), random = ~ 1 | Study.ID, struct="UN", data=df, method=method)
      #     c.	Species loss due to PIE/ENS_PIE: Fragmentation reduces the evenness.
      model[["ENS_PIE"]] <- rma.mv(yi=ES.ENS_pie,V=ES.var.ENS_pie, mods = as.formula(mods.formula), random = ~ 1 | Study.ID, struct="UN", data=df, method=method)
   # 
   # 
      model[["BetaDiv_PA"]] <- rma.mv(yi=repl_part_BS_qF,V=1, mods = as.formula(mods.formula), random = ~ 1 | Study.ID, struct="UN", data=df, method=method)
      model[["BetaDiv_abund"]] <- rma.mv(yi=repl_part_BS_qT,V=1, mods = as.formula(mods.formula), random = ~ 1 | Study.ID, struct="UN", data=df, method=method)
      
   return(model)
}

model_frag <- model_frag_group <- model_gradient <- list()

#--------------------------------------------------------
## Single Responses
# 1. Grand Mean
model_frag[["GrandMean"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~1",method="REML")
model_frag_group[["GrandMean"]] <- analysis_func(df=ES_frag_group_df.complete, mods.formula="~1",method="REML")
model_gradient[["GrandMean"]] <- analysis_func(df=ES_df.complete, mods.formula="~1",method="REML")

# 2. Taxonomic group
#      Do effects vary systematically vary with taxonomic group, i.e. is the overall heterogeneity reduced within groups of the same taxa?
model_frag[["Taxa"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~taxa-1", method="REML")
model_frag_group[["Taxa"]] <- analysis_func(df=ES_frag_group_df.complete, mods.formula="~taxa-1", method="REML")
model_gradient[["Taxa"]] <- analysis_func(df=ES_df.complete, mods.formula="~taxa-1", method="REML")

# 3. Biome
model_frag[["Biome"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~biome-1", method="REML")
model_frag_group[["Biome"]] <- analysis_func(df=ES_frag_group_df.complete, mods.formula="~biome-1", method="REML")
model_gradient[["Biome"]] <- analysis_func(df=ES_df.complete, mods.formula="~biome-1", method="REML")

# 4. Matrix category
model_frag[["Matrix"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~matrix.category-1", method="REML")
model_frag_group[["Matrix"]] <- analysis_func(df=ES_frag_group_df.complete, mods.formula="~matrix.category-1", method="REML")
model_gradient[["Matrix"]] <- analysis_func(df=ES_df.complete, mods.formula="~matrix.category-1", method="REML")

# 5. Time since fragmentation
model_frag[["Time"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~time.since.fragmentation-1", method="REML")
model_frag_group[["Time"]] <- analysis_func(df=ES_frag_group_df.complete, mods.formula="~time.since.fragmentation-1", method="REML")
model_gradient[["Time"]] <- analysis_func(df=ES_df.complete, mods.formula="~time.since.fragmentation-1", method="REML")

# 6. Ratio min/max fragment
model_frag[["SizeRatio"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~ratio.min.max.fragment.size2", method="REML")
model_frag_group[["SizeRatio"]] <- analysis_func(df=ES_frag_group_df.complete, mods.formula="~ratio.min.max.fragment.size2", method="REML")
model_gradient[["SizeRatio"]] <- analysis_func(df=ES_df.complete, mods.formula="~ratio.min.max.fragment.size2", method="REML")

#--------------------------------------------------------
# Additive Effects
# Matrix:Time
model_frag[["Matrix+Time"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~matrix.category+time.since.fragmentation", method="REML")
model_frag_group[["Matrix+Time"]] <- analysis_func(df=ES_frag_group_df.complete,mods.formula="~matrix.category+time.since.fragmentation", method="REML")
model_gradient[["Matrix+Time"]] <- analysis_func(df=ES_df.complete, mods.formula="~matrix.category+time.since.fragmentation", method="REML")

# Matrix+Taxa
model_frag[["Matrix+Taxa"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~matrix.category+taxa", method="REML")
model_frag_group[["Matrix+Taxa"]] <- analysis_func(df=ES_frag_group_df.complete,mods.formula="~matrix.category+taxa", method="REML")
model_gradient[["Matrix+Taxa"]] <- analysis_func(df=ES_df.complete, mods.formula="~matrix.category+taxa", method="REML")

# Matrix+SizeRatio
model_frag[["Matrix+SizeRatio"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~matrix.category+ratio.min.max.fragment.size2", method="REML")
model_frag_group[["Matrix+SizeRatio"]] <- analysis_func(df=ES_frag_group_df.complete,mods.formula="~matrix.category+ratio.min.max.fragment.size2", method="REML")
model_gradient[["Matrix+SizeRatio"]] <- analysis_func(df=ES_df.complete, mods.formula="~matrix.category+ratio.min.max.fragment.size2", method="REML")

# Time+Taxa
model_frag[["Time+Taxa"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~time.since.fragmentation+taxa", method="REML")
model_frag_group[["Time+Taxa"]] <- analysis_func(df=ES_frag_group_df.complete,mods.formula="~time.since.fragmentation+taxa", method="REML")
model_gradient[["Time+Taxa"]] <- analysis_func(df=ES_df.complete, mods.formula="~time.since.fragmentation+taxa", method="REML")
# 

#--------------------------------------------------------
# Interactions
# Matrix:Time
model_frag[["Matrix:Time"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~matrix.category+time.since.fragmentation+matrix.category:time.since.fragmentation", method="REML")
model_frag_group[["Matrix:Time"]] <- analysis_func(df=ES_frag_group_df.complete,mods.formula="~matrix.category+time.since.fragmentation+matrix.category:time.since.fragmentation", method="REML")
model_gradient[["Matrix:Time"]] <- analysis_func(df=ES_df.complete, mods.formula="~matrix.category+time.since.fragmentation+matrix.category:time.since.fragmentation", method="REML")

# Matrix:Taxa
model_frag[["Matrix:Taxa"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~matrix.category+taxa+matrix.category:taxa", method="REML")
model_frag_group[["Matrix:Taxa"]] <- analysis_func(df=ES_frag_group_df.complete,mods.formula="~matrix.category+taxa+matrix.category:taxa", method="REML")
model_gradient[["Matrix:Taxa"]] <- analysis_func(df=ES_df.complete, mods.formula="~matrix.category+taxa+matrix.category:taxa", method="REML")

# Matrix:SizeRatio
model_frag[["Matrix:SizeRatio"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~matrix.category+ratio.min.max.fragment.size2+matrix.category:ratio.min.max.fragment.size2", method="REML")
model_frag_group[["Matrix:SizeRatio"]] <- analysis_func(df=ES_frag_group_df.complete,mods.formula="~matrix.category+ratio.min.max.fragment.size2+matrix.category:ratio.min.max.fragment.size2", method="REML")
model_gradient[["Matrix:SizeRatio"]] <- analysis_func(df=ES_df.complete, mods.formula="~matrix.category+ratio.min.max.fragment.size2+matrix.category:ratio.min.max.fragment.size2", method="REML")

# Time:Taxa
model_frag[["Time:Taxa"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~time.since.fragmentation+taxa+time.since.fragmentation:taxa", method="REML")
model_frag_group[["Time:Taxa"]] <- analysis_func(df=ES_frag_group_df.complete,mods.formula="~time.since.fragmentation+taxa+time.since.fragmentation:taxa", method="REML")
model_gradient[["Time:Taxa"]] <- analysis_func(df=ES_df.complete, mods.formula="~time.since.fragmentation+taxa+time.since.fragmentation:taxa", method="REML")
# 

# sample_design
model_frag[["sample_design"]] <- analysis_func(df=ES_frag_df.complete, mods.formula="~sample_design-1", method="REML")
model_frag_group[["sample_design"]] <- analysis_func(df=ES_frag_group_df.complete,mods.formula="~sample_design-1", method="REML")
model_gradient[["sample_design"]] <- analysis_func(df=ES_df.complete,mods.formula="~sample_design-1", method="REML")


#--------------------------------------------------------
## Model selection for models with interactions
#model_frag[["Matrix:Time"]] <- lapply(model_frag[["Matrix:Time"]], RMAselect) 
## Error in update.rma(currentModel, "~ 1") : 
## The 'mods' argument in 'object' must be a formula for updating to work. 
## 'mods' has to be explicit, not within as.formula() 


#--------------------------------------------------------
### RESTERAMPE
# # 3.	Common vs rare: Fragmentation leads to decrease in coverage standardized species richness and to a change in the species-abundance distribution, e.g. rare species are more affected than common species (increasing ENS_PIE).
# df_long_sub <- subset(df_long, ES %in% c("ES.D0_hat", "ES.ENS_pie"))
# df_long_sub$ES <- factor(df_long_sub$ES)
# model[["Common vs. rare"]] <- rma.mv(yi=value,V=ES.var, mods = as.formula(mods.formula.extended), random = ~ ES | Study.ID, struct="UN", data=df_long_sub, method=method)

