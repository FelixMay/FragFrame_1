load(path2temp %+% "02.1_Data4Analysis_out.Rdata") 
load(path2temp %+% "02.3_DataAnalysis_out.Rdata") # model_frag, model_frag_group, model_gradient

plot.func <- function(model,ylab,title){
   pred.se.list <- lapply(model,function(x) data.frame(b=x$b,se=x$se))
   # pred.se.list <- lapply(model,function(x) lapply(x, function(y) data.frame(b=y$b,se=y$se)))
   #   str(pred.se.list)
   pred.se.df.list <- lapply(pred.se.list, function(x) bind_rows(x,.id="model"))
   pred.se.df <- bind_rows(pred.se.df.list)
#   pred.se.df$levels <- rep(c("frag","frag_group","gradient"),each=4)
   pred.se.df$variable <- c("S", "D0_hat", "N_std", "ENS_pie") # rep(c("S", "D0_hat", "N_std", "ENS_pie"),times=3)
   ylim <- c(floor(10*min(pred.se.df$b-1.96*pred.se.df$se))/10,ceiling(10*max(pred.se.df$b+1.96*pred.se.df$se))/10)
   
   pd <- position_dodge(width=0.4)
#   coul = brewer.pal(length(levels(pred.se.df$levels)), "Accent")    # Classic palette, with n.study colors
   plot1 <- ggplot(pred.se.df) + #, aes(color=levels)) +
      geom_point(position=pd,aes(x=variable, y=b),size=4) +
      geom_errorbar(position=pd,aes(x=variable, ymin=b-1.96*se,ymax=b+1.96*se),width=0.2,size=1.2) +
      geom_hline(yintercept=0,linetype="twodash", size=0.6) +
      scale_x_discrete("",limits=c("S", "D0_hat", "N_std", "ENS_pie")) +
      xlab("") + ylab(ylab) + ylim(ylim) +
      ggtitle(title) +
      # scale_color_manual("", values=coul,
      #                    guide=guide_legend(title="", direction="horizontal",nrow=1)) +
      theme_bw() +
      theme(legend.position="bottom")
   print(plot1)
}

#model <- list(model_frag[["GrandMean"]],model_frag_group[["GrandMean"]],model_gradient[["GrandMean"]])

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_GrandMean.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func(model=model_frag[["GrandMean"]],ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_group_GrandMean.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func(model_frag_group[["GrandMean"]], ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_gradient_GrandMean.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func(model_gradient[["GrandMean"]], ylab="Fishers' z",title="Gradient of fragmentation")
dev.off()

#####################################################################################
plot.func.levels <- function(model,cov,levels,ylab,title){
   pred.se.list <- lapply(model,function(x) data.frame(b=x$b,se=x$se))
   pred.se.df <- bind_rows(pred.se.list,.id="model")  
   pred.se.df$variable <- rep(c("S", "D0_hat", "N_std", "ENS_pie"),each=length(levels))
   pred.se.df$levels <- unlist(lapply(model, function(x) sapply(strsplit(rownames(x$b),cov), function(y) y[2])))
   pred.se.df$levels <- factor(pred.se.df$levels, levels=levels) 
   
   ylim <- c(floor(10*min(pred.se.df$b-1.96*pred.se.df$se))/10,ceiling(10*max(pred.se.df$b+1.96*pred.se.df$se))/10)
   
   pd <- position_dodge(width=0.6)
   coul = brewer.pal(ifelse(length(levels)>3, length(levels), 3), "Accent")    # Classic palette, with n.study colors
   
   plot1 <- ggplot(pred.se.df, aes(color=levels)) +
      geom_point(position=pd,aes(x=variable, y=b),size=4) +
      geom_errorbar(position=pd,aes(x=variable, ymin=b-1.96*se,ymax=b+1.96*se),width=0.2,size=1.2) +
      geom_hline(yintercept=0,linetype="twodash", size=0.6) +
      scale_x_discrete("",limits=c("S", "D0_hat", "N_std", "ENS_pie")) +
      xlab("") + ylab(ylab) + ylim(ylim) +
      ggtitle(title) +
      scale_color_manual("", values=coul,
                         guide=guide_legend(title="", direction="horizontal",nrow=round(length(levels)/3))) +
      theme_bw() +
      theme(legend.position="bottom")
   
   print(plot1)
}

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_Taxa.png", width=30,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model=model_frag[["Taxa"]],cov="taxa", levels=c("invertebrates", "amphibians & reptiles", "birds", "mammals","plants"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_Biome.png", width=30,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model=model_frag[["Biome"]],cov="biome", levels=levels(ES_frag_df.complete$biome), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_Matrix.png", width=30,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model=model_frag[["Matrix"]],cov="matrix.category", levels=c("light filter", "medium filter", "harsh filter"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_Time.png", width=30,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model=model_frag[["Time"]],cov="time.since.fragmentation", levels=c("recent (less than 20 years)", "intermediate (20-100 years)", "long (100+ years)"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

##
png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_group_Taxa.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model_frag_group[["Taxa"]],cov="taxa", levels=c("invertebrates", "amphibians & reptiles", "birds", "mammals","plants"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_group_Biome.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model_frag_group[["Biome"]],cov="biome", levels=levels(ES_frag_df.complete$biome), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_group_Matrix.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model_frag_group[["Matrix"]],cov="matrix.category", levels=c("light filter", "medium filter", "harsh filter"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_group_Time.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model_frag_group[["Time"]],cov="time.since.fragmentation", levels=c("recent (less than 20 years)", "intermediate (20-100 years)", "long (100+ years)"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

##
png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_gradient_Taxa.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model_gradient[["Taxa"]],cov="taxa", levels=c("invertebrates", "amphibians & reptiles", "birds", "mammals","plants"), ylab="Fishers' z",title="Gradient of fragmentation")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_gradient_Biome.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model_gradient[["Biome"]],cov="biome", levels=levels(ES_frag_df.complete$biome), ylab="Fishers' z",title="Gradient of fragmentation")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_gradient_Matrix.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model_gradient[["Matrix"]],cov="matrix.category", levels=c("light filter", "medium filter", "harsh filter"), ylab="Fishers' z",title="Gradient of fragmentation")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_gradient_Time.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model_gradient[["Time"]],cov="time.since.fragmentation", levels=c("recent (less than 20 years)", "intermediate (20-100 years)", "long (100+ years)"), ylab="Fishers' z",title="Gradient of fragmentation")
dev.off()



#####################################################################################
plot.func.continuous <- function(model,newmods,ylab,title){
   pred.se.list <- lapply(model,function(x) predict(x,newmods=newmods))
   #   str(pred.se.list)
   pred.list2df <- function(pred.list){
      pred.df <- data.frame(pred=pred.list$pred,
                            se=pred.list$se,
                            ci.lb=pred.list$ci.lb,
                            ci.ub=pred.list$ci.ub,
                            cr.lb=pred.list$cr.lb,
                            cr.ub=pred.list$cr.ub)
      return(pred.df)
   }
   pred.se.df <- bind_rows(lapply(pred.se.list,pred.list2df),.id="model")
   
   pred.se.df$variable <- factor(rep(c("S", "D0_hat", "N_std", "ENS_pie"),each=length(newmods)),levels=c("S", "D0_hat", "N_std", "ENS_pie"))
   pred.se.df$newmods <- newmods 

   plot1 <- ggplot(pred.se.df,aes(y=pred,x=newmods)) +
      geom_line(size=2) +
      geom_hline(yintercept=0,linetype="twodash", size=0.6) +
      geom_ribbon(aes(ymin=ci.lb,ymax=ci.ub),alpha=0.3,color=NA,fill="grey") +
      xlab("Min/Max Size Ratio") + ylab(ylab) +
      ggtitle(title) +
      facet_grid(variable~.) +
      theme_bw() +
      theme(legend.position="bottom")
   
   print(plot1)
}


png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_SizeRatio.png", width=30,height=10,units="cm",res=200,type = "cairo-png")
plot.func.continuous(model=model_frag[["SizeRatio"]], newmods=seq(min(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T),max(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T),by=0.05), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_group_SizeRatio.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.continuous(model_frag_group[["SizeRatio"]], newmods=seq(min(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T),max(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T),by=0.05), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_gradient_SizeRatio.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.continuous(model_gradient[["SizeRatio"]], newmods=seq(min(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T),max(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T),by=0.05), ylab="Fishers' z",title="Gradient of fragmentation")
dev.off()

#-----------------
png(file=path2temp %+% "ResultsPlots/SensitivityAnalysis/ResultPlot_frag_sample.design.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model_sample.design[["frag"]],cov="sample.design", levels=c("standardized","subsamples_in_frag"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/SensitivityAnalysis/ResultPlot_frag_group_sample.design.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model_sample.design[["frag_group"]],cov="sample.design", levels=c("standardized","pooled","subsamples_in_frag"), ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

png(file=path2temp %+% "ResultsPlots/SensitivityAnalysis/ResultPlot_gradient_sample.design.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func.levels(model_sample.design[["gradient"]],cov="sample.design", levels=c("standardized","pooled","subsamples_in_frag"), ylab="Fishers' z",title="Gradient of fragmentation")
dev.off()

# setRefToMostCommonLevel <- function(f) {
#    f <- as.factor(f)
#    t <- table(f)
#    f.new <- relevel(f,ref=as.integer(which(t>=max(t))[[1]]))
#    return(f.new)
# }
# ES_frag_df.complete$matrix.category <- setRefToMostCommonLevel(ES_frag_df.complete$matrix.category)
# 
# # Matrix:SizeRatio
# model <- model_frag[["Matrix:SizeRatio"]][[1]]
# newdat <- expand.grid(matrix.category=levels(ES_frag_df.complete$matrix.category),ratio.min.max.fragment.size2=seq(min(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T),max(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T),by=0.1))
# mm <- model.matrix(~matrix.category+ratio.min.max.fragment.size2+matrix.category:ratio.min.max.fragment.size2, data=newdat)
# preds <- predict(model,newmods=cbind("harsh filter",seq(min(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T)),by=0.1))
# pred.se.list <- lapply(model,function(x) predict(x,newmods=cbind(mm),intercept=F,addx=T))
# 
# X <- model.matrix(~matrix.category+ratio.min.max.fragment.size2+matrix.category:ratio.min.max.fragment.size2, data=newdat)
# newdat <- expand.grid(matrix.category=levels(setRefToMostCommonLevel(ES_frag_df.complete$matrix.category)),ratio.min.max.fragment.size2=seq(min(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T),max(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T),by=0.1))
# X_new <- model.matrix(~matrix.category+ratio.min.max.fragment.size2+matrix.category:ratio.min.max.fragment.size2, data=newdat)
# 
# preds <- X_new %*% model$b
# fits <- X %*% model$b
# 
# sh2 <- fitted(res)
# 
# # Matrix+SizeRatio
# model <- model_frag[["Matrix+SizeRatio"]][[1]]
# newdat <- expand.grid(matrix.category=levels(ES_frag_df.complete$matrix.category),ratio.min.max.fragment.size2=seq(min(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T),max(ES_frag_df.complete$ratio.min.max.fragment.size2,na.rm=T),by=0.1))
# mm <- model.matrix(~matrix.category+ratio.min.max.fragment.size2, data=newdat)
# preds <- predict(model,newmods=cbind(levels(ES_frag_df.complete$matrix.category),c(0.1,0.2,0.3,0.4,0.5)))
# pred.se.list <- lapply(model,function(x) predict(x,newmods=cbind(mm),intercept=F,addx=T))
# preds <- predict(model,newmods=c(1,0,0),addx=T)
# 
# 

# Matrix:SizeRatio
ES_frag_group_df.complete$matrix.category <- setRefToMostCommonLevel(ES_frag_group_df.complete$matrix.category)

new.SizeRatio <- seq(min(ES_frag_group_df.complete$ratio.min.max.fragment.size2,na.rm=T),max(ES_frag_group_df.complete$ratio.min.max.fragment.size2,na.rm=T),by=0.1)
newdat <- expand.grid(matrix.category=levels(ES_frag_group_df.complete$matrix.category),ratio.min.max.fragment.size2=new.SizeRatio)
mm <- model.matrix(~matrix.category+ratio.min.max.fragment.size2+matrix.category:ratio.min.max.fragment.size2, data=newdat)

pred.se.list <- lapply(model,function(x) predict.rma(x,newmods=mm))
#   str(pred.se.list)
pred.list2df <- function(pred.list){
   pred.df <- data.frame(pred=pred.list$pred,
                         se=pred.list$se,
                         ci.lb=pred.list$ci.lb,
                         ci.ub=pred.list$ci.ub,
                         cr.lb=pred.list$cr.lb,
                         cr.ub=pred.list$cr.ub,
                         matrix.category=newdat$matrix.category,
                         SizeRatio=newdat$ratio.min.max.fragment.size2)
   return(pred.df)
}
pred.se.df <- bind_rows(lapply(pred.se.list,pred.list2df),.id="model")

pred.se.df$variable <- factor(rep(c("S", "D0_hat", "N_std", "ENS_pie"),each=length(new.SizeRatio)*length(levels(ES_frag_group_df.complete$matrix.category))),levels=c("S", "D0_hat", "N_std", "ENS_pie"))

plot1 <- ggplot(pred.se.df,aes(y=pred,x=SizeRatio, color=matrix.category)) +
   geom_line(size=2) +
   geom_hline(yintercept=0,linetype="twodash", size=0.6) +
   geom_ribbon(aes(ymin=ci.lb,ymax=ci.ub,fill=matrix.category),alpha=0.2,color=NA) +
   xlab("Min/Max Size Ratio") + ylab("Log(Response Ratio)") +
   ggtitle("Smallest vs. largest fragment group") +
   facet_grid(variable~.) +
   theme_bw() +
   theme(legend.position="bottom")

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_group_MatrixSizeRatio.png", width=30,height=20,units="cm",res=200,type = "cairo-png")
print(plot1)
dev.off()

#-------------------------------
# Time:Taxa
ES_df.complete$time.since.fragmentation <- setRefToMostCommonLevel(ES_df.complete$time.since.fragmentation)
ES_df.complete$taxa <- setRefToMostCommonLevel(ES_df.complete$taxa)

newdat <- expand.grid(time.since.fragmentation=levels(ES_df.complete$time.since.fragmentation),taxa=levels(ES_df.complete$taxa))
mm <- model.matrix(~time.since.fragmentation+taxa+time.since.fragmentation:taxa, data=newdat)

pred.se.list <- list()

for(i in 1:4){
   model <- model_gradient[["Time:Taxa"]][[i]]
   mm <- mm[,colnames(mm) %in% rownames(model$b)]
   temp <- predict.rma(model,newmods=mm,addx=T)
   pred.se.list[[i]] <- newdat
   pred.se.list[[i]][,c("b","se","ci.lb","ci.ub")] <- cbind(temp$pred, temp$se, temp$ci.lb, temp$ci.ub)
}

pred.se.df <- bind_rows(pred.se.list,.id="model")  
pred.se.df$variable <- rep(c("S", "D0_hat", "N_std", "ENS_pie"),each=nrow(newdat))
pred.se.df$levels <- factor(paste(pred.se.df$time.since.fragmentation,pred.se.df$taxa, sep = "_"))

ylim <- c(floor(10*min(pred.se.df$ci.lb))/10,ceiling(10*max(pred.se.df$ci.ub))/10)

pd <- position_dodge(width=0.6)
coul = colorRampPalette(brewer.pal(5, "Accent"))    # Classic palette, with n.study colors

plot1 <- ggplot(pred.se.df, aes(color=levels)) +
   geom_point(position=pd,aes(x=variable, y=b),size=4) +
   geom_errorbar(position=pd,aes(x=variable, ymin=b-1.96*se,ymax=b+1.96*se),width=0.2,size=1.2) +
   geom_hline(yintercept=0,linetype="twodash", size=0.6) +
   scale_x_discrete("",limits=c("S", "D0_hat", "N_std", "ENS_pie")) +
   xlab("") + ylab("Fisher's z") + ylim(ylim) +
   ggtitle("Gradient of fragmentation") +
 #  guide_legend(title="",direction="horizontal",nrow=round(length(levels)/3)) +
   scale_color_manual("", values=coul(15))+#,
#                      guide=guide_legend(title="", direction="horizontal",nrow=round(length(pred.se.df$levels)/5))) +
   theme_bw() +
   theme(legend.position="bottom")

png(file=path2temp %+% "ResultsPlots/ResultPlot_gradient_TimeTaxa.png", width=30,height=20,units="cm",res=200,type = "cairo-png")
print(plot1)
dev.off()

# dat.bcg$alloc <- factor(dat.bcg$alloc)
# res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, mods = ~ ablat + alloc + ablat:alloc, data=dat.bcg)
# newdat <- expand.grid(ablat=c(10,20,30,40),alloc=levels(dat.bcg$alloc))
# mm <- model.matrix(~ablat + alloc + ablat:alloc, data=newdat)
# preds <- predict(res,newmods=cbind(rep(c(10,20,30,40),times=3),c(rep(0,4),rep(1,4),rep(0,4)),c(rep(0,8),rep(1,4)),rep(c(10,20,30,40),times=3)*c(rep(0,4),rep(1,4),rep(0,4)),rep(c(10,20,30,40),times=3)*c(rep(0,8),rep(1,4))),addx=T)
# predict(res, addx=TRUE)

