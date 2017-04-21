load(path2temp %+% "02.1_Data4Analysis_out.Rdata") 

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
analysis_func <- function(df, covar="intercept.only", method){
   
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
   
   return(model)
}

## Repeat analysis for the 67 studies used in the gradient analysis 
ES_frag_df.complete <- subset(ES_frag_df.complete, n.fragment > 3)
ES_frag_group_df.complete <- subset(ES_frag_group_df.complete, n.fragment > 3)

# 1. Grand Mean
model_frag_sub <- analysis_func(df=ES_frag_df.complete, method="REML")
model_frag_group_sub <- analysis_func(df=ES_frag_group_df.complete, method="REML")

# 2.	Taxonomic group
#      Do effects vary systematically vary with taxonomic group, i.e. is the overall heterogeneity reduced within groups of the same taxa?
modelTaxa_frag_sub <- analysis_func(df=ES_frag_df.complete, covar="taxa", method="REML")
modelTaxa_frag_group_sub <- analysis_func(df=ES_frag_group_df.complete,covar="taxa", method="REML")

plot.func <- function(model, model_sub,ylab,title){
   pred.se.list <- lapply(model,function(x) data.frame(b=x$b,se=x$se))
   pred.se.list_sub <- lapply(model_sub,function(x) data.frame(b=x$b,se=x$se))
   #   str(pred.se.list)
   pred.se.df <- bind_rows(pred.se.list,.id="model")
   pred.se.df_sub <- bind_rows(pred.se.list_sub,.id="model")
   pred.se.df <- rbind(pred.se.df,pred.se.df_sub)
   pred.se.df$data <- rep(c("full", "sub"),each=4)
   pred.se.df$variable <- rep(c("S", "D0_hat", "N_std", "ENS_pie"),times=2)
   ylim <- c(floor(10*min(pred.se.df$b-1.96*pred.se.df$se))/10,ceiling(10*max(pred.se.df$b+1.96*pred.se.df$se))/10)
   
   pd <- position_dodge(width=0.3)
   plot1 <- ggplot(pred.se.df, aes(color=data)) +
      geom_point(position=pd,aes(x=variable, y=b),size=4) +
      geom_errorbar(position=pd,aes(x=variable, ymin=b-1.96*se,ymax=b+1.96*se),width=0.2,size=1.2) +
      geom_hline(yintercept=0,linetype="twodash", size=0.6) +
      scale_x_discrete("",limits=c("S", "D0_hat", "N_std", "ENS_pie")) +
      xlab("") + ylab(ylab) + ylim(ylim) +
      ggtitle(title) +
      theme_habfrag(rel.text.size=1.5)   
   print(plot1)
   
}

plot.func.taxa <- function(model, model_sub,ylab,title){
   pred.se.list <- lapply(model,function(x) data.frame(b=x$b,se=x$se))
   pred.se.list_sub <- lapply(model_sub,function(x) data.frame(b=x$b,se=x$se))
   #   str(pred.se.list)
   pred.se.df <- bind_rows(pred.se.list,.id="model")
   pred.se.df_sub <- bind_rows(pred.se.list_sub,.id="model")
   pred.se.df <- rbind(pred.se.df,pred.se.df_sub)
   pred.se.df$data <- rep(c("full", "sub"),each=4)
   
   pred.se.df$variable <- rep(c("S", "D0_hat", "N_std", "ENS_pie"),each=5)
   pred.se.df$taxa <- c("invertebrates", "amphibians & reptiles", "birds", "mammals","plants") 
   ylim <- c(floor(10*min(pred.se.df$b-1.96*pred.se.df$se))/10,ceiling(10*max(pred.se.df$b+1.96*pred.se.df$se))/10)
   
   pd <- position_dodge(width=0.6)
   
   plot1 <- ggplot(pred.se.df, aes(color=taxa, shape=data)) +
      geom_point(position=pd,aes(x=variable, y=b),size=2) +
      geom_errorbar(position=pd,aes(x=variable, ymin=b-1.96*se,ymax=b+1.96*se),width=0.2,size=1) +
      geom_hline(yintercept=0,linetype="twodash", size=0.6) +
      scale_x_discrete("",limits=c("S", "D0_hat", "N_std", "ENS_pie")) +
      xlab("") + ylab(ylab) + ylim(ylim) +
      ggtitle(title) +
      theme_habfrag(rel.text.size=1.5)   
   
   print(plot1)
}

## Plot results
## uses plot.func from 02.5_VisualizeResults
png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_subset.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func(model=model_frag,model_sub=model_frag_sub, ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_group_subset.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func(model=model_frag_group,model_sub=model_frag_group_sub, ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlotTaxa_frag_subset.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func(model=modelTaxa_frag,model_sub=modelTaxa_frag_sub, ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlotTaxa_frag_group_subset.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.taxa(model=modelTaxa_frag_group,model_sub=modelTaxa_frag_group_sub, ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()
