load(path2temp %+% "02.1_Data4Analysis_out.Rdata") 
load(path2temp %+% "02.3_DataAnalysis_out.Rdata") 

############################################################################
### 2. Study level of consistency among studies
############################################################################

consistency_func <- function(model){
   consistency.tab <- data.frame(model=names(model),
#                                 logLik=NA, deviance=NA, AIC=NA, BIC=NA, AICc=NA,
                                 QE=unlist(lapply(model,function(x) x$QE)), QEp=unlist(lapply(model,function(x) x$QEp)),
                                 QM=unlist(lapply(model,function(x) x$QM)),QMp=unlist(lapply(model,function(x) x$QMp)))
#   consistency.tab[,2:6] <- t(sapply(model,fitstats.rma))
   consistency.tab$R2 <- consistency.tab$QM/(consistency.tab$QM+consistency.tab$QE)
   
   return(consistency.tab)
}

cons_frag <- consistency_func(model_frag)
cons_frag_group <- consistency_func(model_frag_group)
cons_gradient <- consistency_func(model_gradient)

cons_df <- data.frame(data=rep(c("frag","frag_group","gradient"),each=5),bind_rows(list(cons_frag,cons_frag_group,cons_gradient)))

write.csv(cons_df, file=path2temp %+% "model_consistency.csv")

consTaxa_frag <- consistency_func(modelTaxa_frag)
consTaxa_frag_group <- consistency_func(modelTaxa_frag_group)
consTaxa_gradient <- consistency_func(modelTaxa_gradient)

consTaxa_df <- data.frame(data=rep(c("frag","frag_group","gradient"),each=5),bind_rows(list(consTaxa_frag,consTaxa_frag_group,consTaxa_gradient)))

write.csv(consTaxa_df, file=path2temp %+% "modelTaxa_consistency.csv")
############################################################################
### 3. Publication bias
############################################################################
# A funnel plot shows the observed effect sizes or outcomes on the x-axis against some measure of precision of the observed effect sizes or outcomes on the y-axis. Based on Sterne and Egger (2001), the recommended choice for the y-axis is the standard error (in decreasing order) and this is also the default for the funnel() function in the metafor package. In the absence of publication bias and heterogeneity, one would then expect to see the points forming a funnel shape, with the majority of the points falling inside of the pseudo-confidence region with bounds ^θ±1.96SE, where ^θ is the estimated effect or outcome based on the fixed-effects model and SE is the standard error value from the y-axis. With other measures of precision for the y-axis, the expected shape of the funnel can be rather different. (Source: http://www.metafor-project.org/doku.php/plots:funnel_plot_variations)

### Funnel plots with transparent points
funnel.plot <- function(yi,sei,refline,x.label,main.title){
   # function input: yi=vector of observed effect sizes
   #                 sei=vector of the corresponding standard errors
   #                 refline=summary effect
   # function output:funnel plot of effect sizes against precision (i.e. standard error)
   
   alpha <- 0.05
   ylim <- c(0, max(sei) * 1)
   steps <- 5
   
   x.lb.bot <- refline - qnorm(alpha/2, lower.tail = FALSE) * sqrt(ylim[2])
   x.ub.bot <- refline + qnorm(alpha/2, lower.tail = FALSE) * sqrt(ylim[2])
   x.lb.top <- refline - qnorm(alpha/2, lower.tail = FALSE) * sqrt(ylim[1])
   x.ub.top <- refline + qnorm(alpha/2, lower.tail = FALSE) * sqrt(ylim[1])
   
   xlim <- c(min(x.lb.bot, min(yi)), max(x.ub.bot, max(yi)))
   rxlim <- xlim[2] - xlim[1]
   xlim[1] <- xlim[1] - (rxlim * 0.1)
   xlim[2] <- xlim[2] + (rxlim * 0.1)
   
   plot(NA, NA, xlim = xlim, ylim = max(sei) - c(ylim[2], ylim[1]), xlab = x.label, ylab = "Standard Error", 
        xaxt = "n", yaxt = "n", bty = "n",main=main.title)
   
   par.usr <- par("usr")
   rect(par.usr[1], par.usr[3], par.usr[2], par.usr[4], col = "lightgray", border = NA)
   axis(side = 2, at = max(sei) - seq(ylim[2], ylim[1], length = steps), 
        labels = formatC(seq(ylim[2], ylim[1], length = steps), digits = 2, format = "f"),las=1)
   abline(h = max(sei) - seq(ylim[2], ylim[1], length = steps), col = "white")
   
   rylim <- ylim[2] - ylim[1]
   ylim[1] <- max(0, ylim[1] - (rylim * 0.1))
   ylim[2] <- ylim[2] + (rylim * 0.1)
   vi.vals <- seq(ylim[1], ylim[2], length = 100)
   ci.left <- refline - qnorm(alpha/2, lower.tail = FALSE) * sqrt(vi.vals)
   ci.right <- refline + qnorm(alpha/2, lower.tail = FALSE) * sqrt(vi.vals)
   lvi <- length(vi.vals)
   polygon(c(ci.left[lvi:1], ci.right), c(max(sei) - sqrt(vi.vals)[lvi:1], max(sei) - sqrt(vi.vals)), border = NA, col = "white")
   segments(refline, max(sei), refline, max(sei) - ylim[2])
   lines(ci.left, max(sei) - sqrt(vi.vals), lty = "dotted")
   lines(ci.right, max(sei) - sqrt(vi.vals), lty = "dotted")
   box(bty = "l")
   at <- axTicks(side = 1)
   axis(side = 1, at = at, labels = at)
   points(yi, max(sei) - sei,col=rgb(0,0,0,50,maxColorValue=255), pch = 19)
   
}

funnel.plot(residuals(model_gradient[[1]]),sqrt(ES_df.complete$ES.var.S),0,x.label="Residuals",main.title=names(model_gradient)[1])
funnel.plot(residuals(model_gradient[[2]]),sqrt(ES_df.complete$ES.var.D0_hat),0,x.label="Residuals",main.title=names(model_gradient)[2])
funnel.plot(residuals(model_gradient[[3]]),sqrt(ES_df.complete$ES.var.N_std),0,x.label="Residuals",main.title=names(model_gradient)[3])
funnel.plot(residuals(model_gradient[[4]]),sqrt(ES_df.complete$ES.var.ENS_pie),0,x.label="Residuals",main.title=names(model_gradient)[4])
df_long_sub <- subset(ES_df.complete_long, ES %in% c("ES.D0_hat", "ES.ENS_pie"))
funnel.plot(residuals(model_gradient[[5]]),sqrt(df_long_sub$ES.var[as.numeric(names(residuals(model_gradient[[5]])))]),0,x.label="Residuals",main.title=names(model_gradient)[5])# NA effect sizes removed prior to modelling and do not appear in residuals, need to manually remove corresponding variances

# Not sure if this is valuable since we set variances to 1
#funnel.plot(residuals(model_frag[[1]]),rep(1,length(residuals(model_frag[[1]]))),0,x.label="Residuals",main.title=names(model_frag)[1])

### Test of funnel plot asymmetry
## regtest() is not working for rma.mv
asymmetry.test <- function(model,weights){
   residuals.vec <- residuals(model)
   fm <- lm(residuals.vec*sqrt(weights)~sqrt(weights))  
   # test if intercept is significantly different from zero
   if(summary(fm)$coefficients[1,4]<0.05){print("There is evidence for publication bias.")}
   if(summary(fm)$coefficients[1,4]>0.05){print("No evidence for publication bias.")}
   return(summary(fm)$coefficients[1,c(1,4)])
}

asymmetry.test(model_gradient[[1]],1/ES_df.complete$ES.var.S)
asymmetry.test(model_gradient[[2]],1/ES_df.complete$ES.var.D0_hat)
asymmetry.test(model_gradient[[3]],1/ES_df.complete$ES.var.N_std)
asymmetry.test(model_gradient[[4]],1/ES_df.complete$ES.var.ENS_pie)
asymmetry.test(model_gradient[[5]],1/df_long_sub$ES.var[as.numeric(names(residuals(model_gradient[[5]])))])

############################################################################
### 4. Sensitivity analysis
############################################################################
### a plot of Cook's distances against leverage/(1-leverage). 
# influence() is not yet implemented in metafor for rma.mv objects. 
# an alternative approach from http://people.stern.nyu.edu/jsimonof/classes/2301/pdf/diagnost.pdf
influence.func <- function(model){
   h <- hatvalues(model)
   x.seq <- seq(min(rstandard(model$z),max(rstandard(model$z),length.out=100)))
   y.seq <- seq(min(h),max(h),length.out=100)
   D <- (x.seq^2*y.seq)/((length(model$b)+1)*(1-y.seq))
   plot(rstandard(model$z~h,xlab="Leverage",ylab="Standardized Residuals", main="Residuals vs. Leverage"))
                         #contour(x.seq,y.seq,D) # not working
                         
}

influence.func(model_frag[[1]])
