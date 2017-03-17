load(path2temp %+% "Data4Analysis.Rdata") 
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
for(col in c("taxa","country", "continent", "biome")){
   df_long[,col] <- setRefToMostCommonLevel(df_long[,col])
   df[,col] <- setRefToMostCommonLevel(df[,col])
}

############################################################################
### 2. Analysis pipeline
############################################################################
analysis_func <- function(df,df_long){
 
   df <- ES_df.complete #ES_frag_group_df.complete, ES_df.complete
   df_long <- ES_df.complete_long
   
   model <- list()

   if (length(grep("var",df))>0){
      
   # 1.	 Less Individuals
   #     a.	Does the absolute abundance N (standardized by area) changes with fragment size?
      model[["Less individuals"]] <- rma.mv(yi=ES.N_std,V=ES.var.N_std, mods = ~ 1, random = ~ 1 | Study.ID, struct="UN", data=df, method="ML")
      
   # 2.	Species loss
   #     a.	Does the observed species richness (standardized by area) changes with fragment size?
      model[["Species loss a"]] <- rma.mv(yi=ES.S, V=ES.var.S, mods = ~ 1, random = ~ 1 | Study.ID, struct="UN", data=df, method="ML")
      #     b.	Does the extrapolated species richness (standardized by coverage) changes with fragment size?
      model[["Species loss b"]] <- rma.mv(yi=ES.D0_hat, V=ES.var.D0_hat, mods = ~ 1, random = ~ 1 | Study.ID, struct="UN", data=df, method="ML")
      
      # 3.	SAD changes: Fragmentation per se or common vs. rare
      #     a.	Does the observed species richness (standardized by area) changes with fragment size and does PIE remains equal (indicating fragmentation per se) or changes (indicating that common and rare species are differently affected)?
      df_long_sub <- subset(df_long, ES %in% c("ES.S", "ES.ENS_pie"))
      df_long_sub$ES <- factor(df_long_sub$ES)
      model[["SAD changes a"]] <- rma.mv(yi=value,V=ES.var, mods = ~ ES - 1, random = ~ ES | Study.ID, struct="UN", data=df_long_sub, method="ML")
      #     b.	Does the extrapolated species richness (standardized by coverage) changes with fragment size and does PIE remains equal (indicating fragmentation per se) or changes (indicating that common and rare species are differently affected)?
      df_long_sub <- subset(df_long, ES %in% c("ES.D0_hat", "ES.ENS_pie"))
      df_long_sub$ES <- factor(df_long_sub$ES)
      model[["SAD changes b"]] <- rma.mv(yi=value,V=ES.var, mods = ~ ES - 1, random = ~ ES | Study.ID, struct="UN", data=df_long_sub, method="ML")

      # 4.	Taxonomic group
      #     a.	Does effects vary systematically vary with taxonomic group, i.e. is the overall heterogeneity reduced within groups of the same taxa?
      # 5.	Habitat and biome of fragments and matrix
      #     a.	Does effects vary systematically vary with types of habitat and biome of fragments and surrounding matrix, i.e. is the overall heterogeneity reduced within groups of the same habitat/biome type of fragment/matrix or combinations?
   }

   
   if(length(grep("var"),df)==0)
      # 1.	 Less Individuals
      #     a.	Does the absolute abundance N (standardized by area) changes with fragment size?
      model[["Less individuals"]] <- lmer(df[,"ES.N_std"] ~ 1 + (1 | Study.ID), data=df)

      # 2.	Species loss
      #     a.	Does the observed species richness (standardized by area) changes with fragment size?
      model[["Species loss a"]] <- lmer(df[,"ES.S"] ~ 1 + (1 | Study.ID), struct="UN", data=df)
      #     b.	Does the extrapolated species richness (standardized by coverage) changes with fragment size?
      model[["Species loss b"]] <- lmer(df[,"ES.D0_hat"] ~ 1, random = ~ 1 | Study.ID, data=df)

      # 3.	SAD changes: Fragmentation per se or common vs. rare
      #     a.	Does the observed species richness (standardized by area) changes with fragment size and does PIE remains equal (indicating fragmentation per se) or changes (indicating that common and rare species are differently affected)?
      df_long_sub <- subset(df_long, ES %in% c("ES.S", "ES.ENS_pie"))
      df_long_sub$ES <- factor(df_long_sub$ES)
      model[["SAD changes a"]] <- lmer(value ~ ES - 1 + (ES | Study.ID), data=df_long_sub)
      #     b.	Does the extrapolated species richness (standardized by coverage) changes with fragment size and does PIE remains equal (indicating fragmentation per se) or changes (indicating that common and rare species are differently affected)?
      df_long_sub <- subset(df_long, ES %in% c("ES.D0_hat", "ES.ENS_pie"))
      df_long_sub$ES <- factor(df_long_sub$ES)
      model[["SAD changes b"]] <- lmer(value ~ ES - 1 + (ES | Study.ID), data=df_long_sub)
   # 5.	Taxonomic group
   #     a.	Does effects vary systematically vary with taxonomic group, i.e. is the overall heterogeneity reduced within groups of the same taxa?
   # 6.	Habitat and biome of fragments and matrix
   #     a.	Does effects vary systematically vary with types of habitat and biome of fragments and surrounding matrix, i.e. is the overall heterogeneity reduced within groups of the same habitat/biome type of fragment/matrix or combinations?
  
   return(model)
}

save.image(file=path2temp %+% "02.3_DataAnalysis_out.Rdata")

### RESTERAMPE
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

