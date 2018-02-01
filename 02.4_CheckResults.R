load(path2temp %+% "02.1_Data4Analysis_out.Rdata") 
load(path2temp %+% "02.3_DataAnalysis_out.Rdata") 


###########################################################################
### 2. Study level of consistency among studies
############################################################################

consistency_func <- function(model){
   consistency.tab <- data.frame(model=names(model),
                                 n.samples = unlist(lapply(model, function(x) x$k)), # k = number of studies included in the model.
                                 logLik=NA, deviance=NA, AIC=NA, BIC=NA, AICc=NA,
                                 QE=unlist(lapply(model,function(x) x$QE)), QEp=unlist(lapply(model,function(x) x$QEp)),
                                 QM=unlist(lapply(model,function(x) x$QM)),QMp=unlist(lapply(model,function(x) x$QMp)))
   consistency.tab[,3:7] <- t(sapply(model, function(x) fitstats.rma(x, REML=F)))
   consistency.tab$R2 <- consistency.tab$QM/(consistency.tab$QM+consistency.tab$QE)
   
   return(consistency.tab)
}

cons_frag <- bind_rows(lapply(model_frag,consistency_func),.id="Covariate")
cons_frag <- arrange(cons_frag, model, AIC) # sort by response than by AIC

cons_frag_group <- bind_rows(lapply(model_frag_group,consistency_func),.id="Covariate")
cons_frag_group <- arrange(cons_frag_group, model, AIC) # sort by response than by AIC

cons_gradient <- bind_rows(lapply(model_gradient,consistency_func),.id="Covariate")
cons_gradient <- arrange(cons_gradient, model, AIC) # sort by response than by AIC

cons_df <- bind_rows(list(cons_frag,cons_frag_group,cons_gradient), .id="data")
cons_df$data <- factor(cons_df$data, levels=c("1","2","3"), labels=c("frag","frag_group","gradient"))

cons_df[,5:14] <- apply(cons_df[,5:14], 2,function(x) round(x,digits=3))

write.csv(cons_df, file=path2temp %+% "model_consistency.csv")

#---------------------------------------------------------------------------
# Model diagnostics
#---------------------------------------------------------------------------
model.diagnostics.func <- function(model){

   ### identifiability of parameters
   profile(model)
   
   ### Residuals vs fitted
   dat <- join(data.frame(rowID=names(fitted(model)),fit=fitted(model)), data.frame(rowID=names(residuals(model,type="rstandard")),res=rstandard(model)$z),by="rowID")
   plot(dat$fit, dat$res)
   abline(h =0)
   
   ### Pearson model residuals
   plot(residuals(model, type="pearson"))
   abline(h =0)
   
   ### Normality of residuals
   qqnorm(residuals(model))
   qqline(residuals(model))
   print(paste("Testing normality of model residuals:" ))
   print(shapiro.test(residuals(model)))
   
}


#---------------------------------------------------------------------------
# Influence diagnostics
#---------------------------------------------------------------------------
cooks.distance.func <- function(model){
   x <- try(cooks.distance(model))
   if(is.error(x)) return()
   plot(x, ylab="Cook's distance", xlab="Study")
   abline(h=0,lty="dashed")
   segments(x0=1:length(x), y0=0, x1 = 1:length(x), y1 = x)
   return(x)
}

#---------------------------------------------------------------------------
# Publication bias
#---------------------------------------------------------------------------

# A funnel plot shows the observed effect sizes or outcomes on the x-axis against some measure of precision of the observed effect sizes or outcomes on the y-axis. Based on Sterne and Egger (2001), the recommended choice for the y-axis is the standard error (in decreasing order) and this is also the default for the funnel() function in the metafor package. In the absence of publication bias and heterogeneity, one would then expect to see the points forming a funnel shape, with the majority of the points falling inside of the pseudo-confidence region with bounds ^θ±1.96SE, where ^θ is the estimated effect or outcome based on the fixed-effects model and SE is the standard error value from the y-axis. With other measures of precision for the y-axis, the expected shape of the funnel can be rather different. (Source: http://www.metafor-project.org/doku.php/plots:funnel_plot_variations)

pub.bias.func <- function(model){
   ### funnel plot
   ### using pre-implemented funnel function
   funnel(residuals(model),vi=model$vi, yaxis="seinv", level=c(90, 95, 99), ylab="Precision (1/SE)",
          xlab="Model residuals", shade=c("white", "gray", "darkgray"), refline=0)
   
   ### asymmetry test
   ### Test of funnel plot asymmetry using Egger's regression (eqn 36 in Nakagawa & Santos 2012 EvolEcol)
   ## regtest() is not working for rma.mv
   residuals.vec <- residuals(model)
   weights <- 1/model$vi
   fm <- lm(residuals.vec*sqrt(weights)~sqrt(weights))
   # test if intercept is significantly different from zero
   if(summary(fm)$coefficients[1,4]<0.05){print("There is evidence for publication bias.")}
   if(summary(fm)$coefficients[1,4]>0.05){print("No evidence for publication bias.")}
   print(summary(fm)$coefficients[1,c(1,4)])
}

#---------------------------------------------------------------------------
# Loop over all models and save prints into a pdf
#---------------------------------------------------------------------------
for(i in names(model_frag)){
   for(j in names(model_frag[[i]])){
      
      print(paste(i, j))
      
      pdf(path2temp %+% "CheckResults/ModelDiagnostics_frag_" %+% i %+% "_" %+% j %+% ".pdf")
      model.diagnostics.func(model_frag[[i]][[j]])
      cooks.distance.func(model_frag[[i]][[j]])
      pub.bias.func(model_frag[[i]][[j]])
      dev.off()
      
      pdf(path2temp %+% "CheckResults/ModelDiagnostics_frag_group_" %+% i %+% "_" %+% j %+% ".pdf")
      model.diagnostics.func(model_frag_group[[i]][[j]])
      cooks.distance.func(model_frag_group[[i]][[j]])
      pub.bias.func(model_frag_group[[i]][[j]])
      dev.off()
      
      pdf(path2temp %+% "CheckResults/ModelDiagnostics_gradient_" %+% i %+% "_" %+% j %+% ".pdf")
      model.diagnostics.func(model_gradient[[i]][[j]])
      cooks.distance.func(model_gradient[[i]][[j]])
      pub.bias.func(model_gradient[[i]][[j]])
      dev.off()
      
   }
}
