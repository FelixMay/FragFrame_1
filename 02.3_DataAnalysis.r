load(path2temp %+% "Data4Analysis.Rdata") # ES_df.complete, df.complete

############################################################################
### 1. set reference levels to the most common levels
############################################################################
setRefToMostCommonLevel <- function(f) {
   f <- as.factor(f)
   t <- table(f)
   f.new <- relevel(f,ref=as.integer(which(t>=max(t))[[1]]))
   return(f.new)
}
df.complete[,c("taxa","country", "continent", "biome")] <- sapply(df.complete[,c("taxa","country", "continent", "biome")],setRefToMostCommonLevel)
ES_df.complete[,c("taxa","country", "continent", "biome")] <- sapply(ES_df.complete[,c("taxa","country", "continent", "biome")],setRefToMostCommonLevel)

############################################################################
### 2. univariate model
############################################################################
model.uv <- list()
model.uv[["N"]] <- rma.mv(z.N, z.var.N, mods = ~ 1, random = ~ 1 | Study.ID, struct="UN", data=ES_df.complete, method="ML")
model.uv[["N + taxa"]] <- rma.mv(z.N, z.var.N, mods = ~ taxa - 1, random = ~ 1 | Study.ID, struct="UN", data=ES_df.complete, method="ML")
model.uv[["S_cov"]] <- rma.mv(z.S_cov, z.var.S_cov, mods = ~ 1, random = ~ 1 | Study.ID, struct="UN", data=ES_df.complete, method="ML")
model.uv[["S_cov + taxa"]] <- rma.mv(z.S_cov, z.var.S_cov, mods = ~ taxa - 1, random = ~ 1 | Study.ID, struct="UN", data=ES_df.complete, method="ML")
model.uv[["PIE"]] <- rma.mv(z.PIE, z.var.PIE, mods = ~ 1, random = ~ 1 | Study.ID, struct="UN", data=ES_df.complete, method="ML")
model.uv[["PIE + taxa"]] <- rma.mv(z.PIE, z.var.PIE, mods = ~ taxa - 1, random = ~ 1 | Study.ID, struct="UN", data=ES_df.complete, method="ML")

fit.tab.uv <- data.frame(model=names(model.uv),
                      logLik=NA, deviance=NA, AIC=NA, BIC=NA, AICc=NA, 
                      QE=unlist(lapply(model.uv,function(x) x$QE)), QEp=unlist(lapply(model.uv,function(x) x$QEp)),
                      QM=unlist(lapply(model.uv,function(x) x$QM)),QMp=unlist(lapply(model.uv,function(x) x$QMp)))
fit.tab.uv[,2:6] <- t(sapply(model.uv,fitstats.rma))
fit.tab.uv$R2 <- fit.tab.uv$QM/(fit.tab.uv$QM+fit.tab.uv$QE)

pred.func <- function(model, newdat, tau2.levels=NULL){
   if(nrow(newdat)==1){
      newdat[,c("z","z.se","z.ci.lb","z.ci.ub")] <- cbind(model$b,model$se,model$b - (1.96*model$se),model$b + (1.96*model$se))
      return(newdat)
   }
   if (nrow(newdat)>1){
      mm <- model.matrix(as.formula(paste("~ - 1 +", names(newdat),collapse="+")),data=newdat)
      preds <- predict.rma(model, newmods = mm, tau2.levels = newdat[,tau2.levels])
      newdat[,c("z","z.se","z.ci.lb","z.ci.ub")] <- cbind(preds$pred,preds$se,preds$ci.lb,preds$ci.ub)   
      return(newdat)
   }
}

preds.uv <- list()
preds.uv[["N"]] <- pred.func(model=model.uv[["N"]],newdat=data.frame(1))
preds.uv[["S_cov"]] <- pred.func(model=model.uv[["S_cov"]],newdat=data.frame(1))
preds.uv[["PIE"]] <- pred.func(model=model.uv[["PIE"]],newdat=data.frame(1))
preds.uv[["N + taxa"]] <- pred.func(model=model.uv[["N + taxa"]],newdat=expand.grid(taxa=levels(ES_df.complete$taxa)))
preds.uv[["S_cov + taxa"]] <- pred.func(model=model.uv[["S_cov + taxa"]],newdat=expand.grid(taxa=levels(ES_df.complete$taxa)))
preds.uv[["PIE + taxa"]] <- pred.func(model=model.uv[["PIE + taxa"]],newdat=expand.grid(taxa=levels(ES_df.complete$taxa)))

############################################################################
### 2. multivariate model 
### cf. http://www.metafor-project.org/doku.php/analyses:vanhouwelingen2002#bivariate_approach
############################################################################
model.mv <- list()
model.mv[["BD"]] <- rma.mv(z, z.var, mods = ~ BD - 1, random = ~ BD | Study.ID, struct="UN", data=df.complete, method="ML")
model.mv[["BD + taxa"]] <- rma.mv(z, z.var, mods = ~ BD + taxa -1, random = ~ BD | Study.ID, struct="UN", data=df.complete, method="ML")

fit.tab.mv <- data.frame(model=names(model.mv),
                         logLik=NA, deviance=NA, AIC=NA, BIC=NA, AICc=NA, 
                         QE=unlist(lapply(model.mv,function(x) x$QE)), QEp=unlist(lapply(model.mv,function(x) x$QEp)),
                         QM=unlist(lapply(model.mv,function(x) x$QM)),QMp=unlist(lapply(model.mv,function(x) x$QMp)))
fit.tab.mv[,2:6] <- t(sapply(model.mv,fitstats.rma))
fit.tab.mv$R2 <- fit.tab.mv$QM/(fit.tab.mv$QM+fit.tab.mv$QE)

preds.mv <- list()
preds.mv[["BD"]] <- pred.func(model=model.mv[["BD"]],newdat=expand.grid(BD=levels(df.complete$BD)),tau2.levels = "BD")
preds.mv[["BD + taxa"]] <- pred.func(model=model.mv[["BD + taxa"]],newdat=expand.grid(BD=levels(df.complete$BD), taxa=levels(factor(df.complete$taxa))),tau2.levels = "BD")

save.image(file=path2temp %+% "02.3_DataAnalysis_out.Rdata")
