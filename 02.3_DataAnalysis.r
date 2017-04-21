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
analysis_func <- function(df,df_long, covar="intercept.only", method){
 
   # df <- ES_df.complete #ES_frag_group_df.complete, ES_df.complete
   # df_long <- ES_df.complete_long
   
   ### set reference levels
   for(col in c("taxa","country", "continent", "biome", "fragment.biome","matrix.biome", "fragment.veg", "matrix.veg")){
      df_long[,col] <- setRefToMostCommonLevel(df_long[,col])
      df[,col] <- setRefToMostCommonLevel(df[,col])
   }
   
   model <- list()

   if(covar=="intercept.only"){
      mods.formula <- "~1"
      mods.formula.extended <- "~ES-1"
   }
   if(covar %in% c("taxa", "country", "continent", "biome", "fragment.biome","matrix.biome", "fragment.veg", "matrix.veg" )){
      mods.formula <- "~" %+% covar %+% "-1"
      mods.formula.extended <- "~ES:" %+% covar %+% "-1"
      
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
      df_long <- data.frame(df_long,
                       ES.var = rep(1,nrow(df)))
      
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
      
      # 3.	Common vs rare: Fragmentation leads to decrease in coverage standardized species richness and to a change in the species-abundance distribution, e.g. rare species are more affected than common species (increasing ENS_PIE).
      df_long_sub <- subset(df_long, ES %in% c("ES.D0_hat", "ES.ENS_pie"))
      df_long_sub$ES <- factor(df_long_sub$ES)
      model[["Common vs. rare"]] <- rma.mv(yi=value,V=ES.var, mods = as.formula(mods.formula.extended), random = ~ 1 | Study.ID, struct="UN", data=df_long_sub, method=method)
   return(model)
}

# 1. Grand Mean
model_frag <- analysis_func(df=ES_frag_df.complete,df_long=ES_frag_df.complete_long, method="REML")
model_frag_group <- analysis_func(df=ES_frag_group_df.complete,df_long=ES_frag_group_df.complete_long, method="REML")
model_gradient <- analysis_func(df=ES_df.complete,df_long=ES_df.complete_long, method="REML")

# 2.	Taxonomic group
#      Do effects vary systematically vary with taxonomic group, i.e. is the overall heterogeneity reduced within groups of the same taxa?
modelTaxa_frag <- analysis_func(df=ES_frag_df.complete,df_long=ES_frag_df.complete_long, covar="taxa", method="REML")
modelTaxa_frag_group <- analysis_func(df=ES_frag_group_df.complete,df_long=ES_frag_group_df.complete_long, covar="taxa", method="REML")
modelTaxa_gradient <- analysis_func(df=ES_df.complete,df_long=ES_df.complete_long, covar="taxa", method="REML")

# 3.	Habitat and biome of fragments and matrix
#     Do effects vary systematically vary with types of habitat and biome of fragments and surrounding matrix, i.e. is the overall heterogeneity reduced within groups of the same habitat/biome type of fragment/matrix or combinations?
### TO DO: Come up with useful categories of habitat and biome of fragments and matrix


### RESTERAMPE
# # ############################################################################
# # 1.	Species loss: Fragmentation leads to decrease of observed species richness, i.e. naïve species richness without considering methodological differences.  
# modelTaxa[["Naive species loss"]] <- rma.mv(yi=ES.S, V=ES.var.S, mods = ~ taxa - 1, random = ~ 1 | Study.ID, struct="UN", data=df, method="REML")
# 
# # 2.	Fragmentation per se:
# #    a.	Species loss due to less coverage standardized species richness: Fragmentation leads to decrease in coverage standardized species richness.
# modelTaxa[["Coverage standardized species loss"]] <- rma.mv(yi=ES.D0_hat, V=ES.var.D0_hat, mods = ~ taxa - 1, random = ~ 1 | Study.ID, struct="UN", data=df, method="REML")
# 
# # b.	Species loss due to less Individuals: Fragmentation decreases the number of individuals.
# modelTaxa[["Less individuals"]] <- rma.mv(yi=ES.N_std,V=ES.var.N_std, mods = ~ taxa - 1, random = ~ 1 | Study.ID, struct="UN", data=df, method="REML")
# 
# # c.	Species loss due to PIE/ENS_PIE: Fragmentation reduces the evenness.
# modelTaxa[["Lower evenness"]] <- rma.mv(yi=ES.ENS_pie,V=ES.var.ENS_pie, mods = ~ taxa - 1, random = ~ 1 | Study.ID, struct="UN", data=df, method="REML")
# 
# # 3.	Common vs rare: Fragmentation leads to decrease in coverage standardized species richness and to a change in the species-abundance distribution, e.g. rare species are more affected than common species (increasing ENS_PIE).
# df_long_sub <- subset(df_long, ES %in% c("ES.D0_hat", "ES.ENS_pie"))
# df_long_sub$ES <- factor(df_long_sub$ES)
# modelTaxa[["Common vs. rare"]] <- rma.mv(yi=value,V=ES.var, mods = ~ ES:taxa - 1, random = ~ ES | Study.ID, struct="UN", data=df_long_sub, method="REML")

# ### unweighted meta-analysis
# if(!any(grepl("var",names(df)))){
#    dfvar <- data.frame(df,var = rep(1,nrow(df)))
# # 1.	Species loss: Fragmentation leads to decrease of observed species richness, i.e. naïve species richness without considering methodological differences.  
#    model[["naive species loss"]] <- lmer(df[,"ES.S"] ~ 1 + (1 | Study.ID), data=df, REML=F)
# # 2.	Fragmentation per se:
#    #    a.	Species loss due to less coverage standardized species richness: Fragmentation leads to decrease in coverage standardized species richness.
#    model[["covstd species loss"]] <- lmer(df[,"ES.D0_hat"] ~ 1 + (1 | Study.ID), data=df, REML=F)
#       # b.	Species loss due to less Individuals: Fragmentation decreases the number of individuals.
#    model[["Less individuals"]] <- lmer(df[,"ES.N_std"] ~ 1 + (1 | Study.ID), data=df, REML=F)
#       # c.	Species loss due to PIE/ENS_PIE: Fragmentation reduces the evenness.
#    model[["Lower evenness"]] <- lmer(df[,"ES.ENS_pie"] ~ 1 + (1 | Study.ID), data=df, REML=F)
#    
#    # 3.	Common vs rare: Fragmentation leads to decrease in coverage standardized species richness and to a change in the species-abundance distribution, e.g. rare species are more affected than common species (increasing ENS_PIE).
#    df_long_sub <- subset(df_long, ES %in% c("ES.D0_hat", "ES.ENS_pie"))
#    df_long_sub$ES <- factor(df_long_sub$ES)
#    model[["Common vs rare"]] <- lmer(value ~ ES - 1 + (ES | Study.ID), data=df_long_sub, REML=F)
#    
#    # 4.	Taxonomic groups: Strength of mechanisms may differ between taxonomic groups. E.g.,….
#    # 5.	Habitat and biome type of the fragment and the surrounding matrix: Strength of mechanisms may differ. Effects might be stronger if natural fragments are surrounded by an intensively managed matrix. 
# }

# ############################################################################
# ### 2. univariate model
# ############################################################################
# model.uv <- list()
# for( BD in BDmetrics){
#    model.uv[[BD]] <- rma.mv(ES_df.complete[,paste("z.",BD,sep="")], ES_df.complete[,paste("z.var.",BD,sep="")], mods = ~ 1, random = ~ 1 | Study.ID, struct="UN", data=ES_df.complete, method="ML")
#    model.uv[[paste(BD,"taxa",sep="+")]] <- rma.mv(ES_df.complete[,paste("z.",BD,sep="")], ES_df.complete[,paste("z.var.",BD,sep="")], mods = ~ taxa-1, random = ~ 1 | Study.ID, struct="UN", data=ES_df.complete, method="ML")
# }
# 
# fit.tab.uv <- data.frame(model=names(model.uv),
#                       logLik=NA, deviance=NA, AIC=NA, BIC=NA, AICc=NA, 
#                       QE=unlist(lapply(model.uv,function(x) x$QE)), QEp=unlist(lapply(model.uv,function(x) x$QEp)),
#                       QM=unlist(lapply(model.uv,function(x) x$QM)),QMp=unlist(lapply(model.uv,function(x) x$QMp)))
# fit.tab.uv[,2:6] <- t(sapply(model.uv,fitstats.rma))
# fit.tab.uv$R2 <- fit.tab.uv$QM/(fit.tab.uv$QM+fit.tab.uv$QE)
# 
# pred.func <- function(model, newdat, tau2.levels=NULL){
#    if(nrow(newdat)==1){
#       newdat[,c("z","z.se","z.ci.lb","z.ci.ub")] <- cbind(model$b,model$se,model$b - (1.96*model$se),model$b + (1.96*model$se))
#       return(newdat)
#    }
#    if (nrow(newdat)>1){
#       mm <- model.matrix(as.formula(paste("~ - 1 +", names(newdat),collapse="+")),data=newdat)
#       preds <- predict.rma(model, newmods = mm, tau2.levels = newdat[,tau2.levels])
#       newdat[,c("z","z.se","z.ci.lb","z.ci.ub")] <- cbind(preds$pred,preds$se,preds$ci.lb,preds$ci.ub)   
#       return(newdat)
#    }
# }
# 
# preds.uv <- list()
# for( BD in BDmetrics){
#    preds.uv[[BD]] <- pred.func(model=model.uv[[BD]],newdat=data.frame(1))
#    preds.uv[[paste(BD,"taxa",sep="+")]] <- pred.func(model=model.uv[[paste(BD,"taxa",sep="+")]],newdat=expand.grid(taxa=levels(ES_df.complete$taxa)))
# }
# 
# ############################################################################
# ### 2. multivariate model 
# ### cf. http://www.metafor-project.org/doku.php/analyses:vanhouwelingen2002#bivariate_approach
# ############################################################################
# model.mv <- list()
# model.mv[["BD"]] <- rma.mv(z, z.var, mods = ~ BD - 1, random = ~ BD | Study.ID, struct="UN", data=df.complete, method="ML")
# model.mv[["BD + taxa"]] <- rma.mv(z, z.var, mods = ~ BD + taxa -1, random = ~ BD | Study.ID, struct="UN", data=df.complete, method="ML")
# 
# fit.tab.mv <- data.frame(model=names(model.mv),
#                          logLik=NA, deviance=NA, AIC=NA, BIC=NA, AICc=NA, 
#                          QE=unlist(lapply(model.mv,function(x) x$QE)), QEp=unlist(lapply(model.mv,function(x) x$QEp)),
#                          QM=unlist(lapply(model.mv,function(x) x$QM)),QMp=unlist(lapply(model.mv,function(x) x$QMp)))
# fit.tab.mv[,2:6] <- t(sapply(model.mv,fitstats.rma))
# fit.tab.mv$R2 <- fit.tab.mv$QM/(fit.tab.mv$QM+fit.tab.mv$QE)
# 
# preds.mv <- list()
# preds.mv[["BD"]] <- pred.func(model=model.mv[["BD"]],newdat=expand.grid(BD=levels(df.complete$BD)),tau2.levels = "BD")
# preds.mv[["BD + taxa"]] <- pred.func(model=model.mv[["BD + taxa"]],newdat=expand.grid(BD=levels(df.complete$BD), taxa=levels(factor(df.complete$taxa))),tau2.levels = "BD")

