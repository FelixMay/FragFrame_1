load(path2temp %+% "02.1_Data4Analysis_out.Rdata") 
load(path2temp %+% "02.3_DataAnalysis_out.Rdata") # model_frag, model_frag_group, model_gradient

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
   for(col in c("taxa","biome","matrix.category","time.since.fragmentation")){
      df[,col] <- setRefToMostCommonLevel(df[,col])
   }
   
   model <- list()
   
   if(covar=="intercept.only"){
      mods.formula <- "~1"
   }
   if(covar %in% c("taxa","biome","matrix.category","time.since.fragmentation")){
      mods.formula <- "~" %+% covar %+% "-1"
   }
   if(covar == "ratio.min.max.fragment.size2"){
      mods.formula <- "~" %+% covar
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

## Repeat analysis for the 58 studies used in the gradient analysis 
ES_frag_df.complete <- subset(ES_frag_df.complete, n.fragment > 3)
ES_frag_group_df.complete <- subset(ES_frag_group_df.complete, n.fragment > 3)

model_frag_sub <- model_frag_group_sub<- list()

# 1. Grand Mean
model_frag_sub[["GrandMean"]] <- analysis_func(df=ES_frag_df.complete, method="REML")
model_frag_group_sub[["GrandMean"]] <- analysis_func(df=ES_frag_group_df.complete, method="REML")

# 2. Taxonomic group
#      Do effects vary systematically vary with taxonomic group, i.e. is the overall heterogeneity reduced within groups of the same taxa?
model_frag_sub[["Taxa"]] <- analysis_func(df=ES_frag_df.complete, covar="taxa", method="REML")
model_frag_group_sub[["Taxa"]] <- analysis_func(df=ES_frag_group_df.complete,covar="taxa", method="REML")

# 3. Biome
model_frag_sub[["Biome"]] <- analysis_func(df=ES_frag_df.complete, covar="biome", method="REML")
model_frag_group_sub[["Biome"]] <- analysis_func(df=ES_frag_group_df.complete,covar="biome", method="REML")

# 4. Matrix category
model_frag_sub[["Matrix"]] <- analysis_func(df=ES_frag_df.complete, covar="matrix.category", method="REML")
model_frag_group_sub[["Matrix"]] <- analysis_func(df=ES_frag_group_df.complete,covar="matrix.category", method="REML")

# 5. Time since fragmentation
model_frag_sub[["Time"]] <- analysis_func(df=ES_frag_df.complete, covar="time.since.fragmentation", method="REML")
model_frag_group_sub[["Time"]] <- analysis_func(df=ES_frag_group_df.complete,covar="time.since.fragmentation", method="REML")

# 6. Ratio min/max fragment
model_frag_sub[["SizeRatio"]] <- analysis_func(df=ES_frag_df.complete, covar="ratio.min.max.fragment.size2", method="REML")
model_frag_group_sub[["SizeRatio"]] <- analysis_func(df=ES_frag_group_df.complete,covar="ratio.min.max.fragment.size2", method="REML")

plot.func <- function(model,model_sub,ylab,title){
   pred.se.list <- lapply(model,function(x) data.frame(b=x$b,se=x$se))
   pred.se.list_sub <- lapply(model_sub,function(x) data.frame(b=x$b,se=x$se))
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
      theme(axis.text = element_text(size = rel(0.8)), 
            axis.ticks = element_line(colour = "black"), 
            axis.title = element_text(size = rel(0.8)), 
            panel.background = element_rect(fill = "white", colour = NA), 
            panel.border = element_rect(fill = NA, colour = "grey50"), 
            panel.grid.major = element_line(colour = "grey90", size = 0.2), 
            panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
            strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2),
            strip.text = element_text(size=rel(0.8)))
   print(plot1)

}

png(file=path2temp %+% "ResultsPlots/SensitivityAnalysis/ResultPlot_frag_GrandMean_subset.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func(model=model_frag[["GrandMean"]],model_sub=model_frag_sub[["GrandMean"]],ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/SensitivityAnalysis/ResultPlot_frag_group_GrandMean_subset.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func(model=model_frag_group[["GrandMean"]], model_sub=model_frag_group_sub[["GrandMean"]], ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

plot.func.levels <- function(model,model_sub,cov,levels,ylab,title){
   pred.se.list <- lapply(model,function(x) data.frame(b=x$b,se=x$se))
   pred.se.list_sub <- lapply(model_sub,function(x) data.frame(b=x$b,se=x$se))
   pred.se.df <- bind_rows(pred.se.list,.id="model")
   pred.se.df_sub <- bind_rows(pred.se.list_sub,.id="model")
   pred.se.df <- rbind(pred.se.df,pred.se.df_sub)
   pred.se.df$data <- rep(c("full", "sub"),each=4)
   pred.se.df$variable <- rep(c("S", "D0_hat", "N_std", "ENS_pie"),each=length(levels))
   pred.se.df$levels <- c(unlist(lapply(model, function(x) sapply(strsplit(rownames(x$b),cov), function(y) y[2]))),
                          unlist(lapply(model_sub, function(x) sapply(strsplit(rownames(x$b),cov), function(y) y[2]))))
   pred.se.df$levels <- factor(pred.se.df$levels, levels=levels) 
   ylim <- c(floor(10*min(pred.se.df$b-1.96*pred.se.df$se))/10,ceiling(10*max(pred.se.df$b+1.96*pred.se.df$se))/10)
   
   pd <- position_dodge(width=0.6)
   coul = brewer.pal(length(levels), "Accent")    # Classic palette, with n.study colors
   
   plot1 <- ggplot(pred.se.df, aes(color=levels, shape=data)) +
      geom_point(position=pd,aes(x=variable, y=b),size=2) +
      geom_errorbar(position=pd,aes(x=variable, ymin=b-1.96*se,ymax=b+1.96*se),width=0.2,size=1) +
      geom_hline(yintercept=0,linetype="twodash", size=0.6) +
      scale_x_discrete("",limits=c("S", "D0_hat", "N_std", "ENS_pie")) +
      xlab("") + ylab(ylab) + ylim(ylim) +
      ggtitle(title) +
      scale_color_manual("", values=coul,
                         guide=guide_legend(title="", direction="horizontal",nrow=round(length(levels)/3))) +
      theme(axis.text = element_text(size = rel(0.8)), 
            axis.ticks = element_line(colour = "black"), 
            axis.title = element_text(size = rel(0.8)), 
            panel.background = element_rect(fill = "white", colour = NA), 
            panel.border = element_rect(fill = NA, colour = "grey50"), 
            panel.grid.major = element_line(colour = "grey90", size = 0.2), 
            panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
            strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2),
            strip.text = element_text(size=rel(0.8)),
            legend.key = element_rect(fill = "white", colour = NA),
            legend.position="bottom")
   
   print(plot1)
}

png(file=path2temp %+% "ResultsPlots/SensitivityAnalysis/ResultPlot_frag_Taxa_subset.png", width=30,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model=model_frag[["Taxa"]], model_sub=model_frag_sub[["Taxa"]],cov="taxa",levels=c("invertebrates", "amphibians & reptiles", "birds", "mammals","plants"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/SensitivityAnalysis/ResultPlot_frag_Biome_subset.png", width=30,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model=model_frag[["Biome"]], model_sub=model_frag_sub[["Biome"]],cov="biome", levels=levels(ES_frag_df.complete$biome), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/SensitivityAnalysis/ResultPlot_frag_Matrix_subset.png", width=30,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model=model_frag[["Matrix"]], model_sub=model_frag_sub[["Matrix"]],cov="matrix.category",levels=c("light filter", "medium filter", "harsh filter"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/SensitivityAnalysis/ResultPlot_frag_Time_subset.png", width=30,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model=model_frag[["Time"]], model_sub=model_frag_sub[["Time"]],cov="time.since.fragmentation",levels=c("recent (less than 20 years)", "intermediate (20-100 years)", "long (100+ years)"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

##
png(file=path2temp %+% "ResultsPlots/SensitivityAnalysis/ResultPlot_frag_group_Taxa_subset.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model=model_frag_group[["Taxa"]], model_sub=model_frag_group_sub[["Taxa"]],cov="taxa",levels=c("invertebrates", "amphibians & reptiles", "birds", "mammals","plants"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_group_Biome_subset.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model=model_frag_group[["Biome"]], model_sub=model_frag_group_sub[["Biome"]],cov="biome",levels=levels(ES_frag_df.complete$biome), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

png(file=path2temp %+% "ResultsPlots/SensitivityAnalysis/ResultPlot_frag_group_Matrix_subset.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model=model_frag_group[["Matrix"]], model_sub=model_frag_group_sub[["Matrix"]],cov="matrix.category",levels=c("light filter", "medium filter", "harsh filter"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

png(file=path2temp %+% "ResultsPlots/SensitivityAnalysis/ResultPlot_frag_group_Time_subset.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model=model_frag_group[["Time"]], model_sub=model_frag_group_sub[["Time"]],cov="time.since.fragmentation", levels=c("recent (less than 20 years)", "intermediate (20-100 years)", "long (100+ years)"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

# TO DO: Sensitivity analysis for Size ratio
# plot.func.continuous <- function(model,newmods,ylab,title){
#    pred.se.list <- lapply(model,function(x) predict(x,newmods=newmods))
#    #   str(pred.se.list)
#    pred.list2df <- function(pred.list){
#       pred.df <- data.frame(pred=pred.list$pred,
#                             se=pred.list$se,
#                             ci.lb=pred.list$ci.lb,
#                             ci.ub=pred.list$ci.ub,
#                             cr.lb=pred.list$cr.lb,
#                             cr.ub=pred.list$cr.ub)
#       return(pred.df)
#    }
#    pred.se.df <- bind_rows(lapply(pred.se.list,pred.list2df),.id="model")
#    
#    pred.se.df$variable <- factor(rep(c("S", "D0_hat", "N_std", "ENS_pie"),each=length(newmods)),levels=c("S", "D0_hat", "N_std", "ENS_pie"))
#    pred.se.df$newmods <- newmods 
#    
#    plot1 <- ggplot(pred.se.df,aes(y=pred,x=newmods)) +
#       geom_line(size=2) +
#       geom_hline(yintercept=0,linetype="twodash", size=0.6) +
#       geom_ribbon(aes(ymin=ci.lb,ymax=ci.ub),alpha=0.3,color=NA,fill="grey") +
#       xlab("Min/Max Size Ratio") + ylab(ylab) +
#       ggtitle(title) +
#       facet_grid(variable~.) +
#       theme(axis.text = element_text(size = rel(0.8)), 
#             axis.ticks = element_line(colour = "black"), 
#             axis.title = element_text(size = rel(0.8)), 
#             panel.background = element_rect(fill = "white", colour = NA), 
#             panel.border = element_rect(fill = NA, colour = "grey50"), 
#             panel.grid.major = element_line(colour = "grey90", size = 0.2), 
#             panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
#             strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2),
#             strip.text = element_text(size=rel(0.8)),
#             legend.key = element_rect(fill = "white", colour = NA),
#             legend.position="bottom")
#    
#    print(plot1)
# }

