load(path2temp %+% "Data4Analysis.Rdata") # ES_df.complete, df.complete
BDmetrics <- c("N", "S", "S_cov", "PIE")

############################################################################
### 1. Plot effect sizes in a 2-dim space
############################################################################
pdf(file=path2temp %+% "2D_Scatterplot_of_z.pdf")
ggplot(data=ES_df.complete, aes(x=z.N,y=z.S,col=taxa)) + 
   geom_point(alpha=.5,size=3) +
   geom_errorbar(aes(ymin=z.S - (1.96*sqrt(z.var.S)),ymax=z.S + (1.96*sqrt(z.var.S)))) +
   geom_errorbarh(aes(xmin=z.N - (1.96*sqrt(z.var.N)),xmax=z.N + (1.96*sqrt(z.var.N)))) +
   geom_abline(intercept=0,slope=1,lty="dotted")

ggplot(data=ES_df.complete, aes(x=z.N,y=z.S_cov,col=taxa)) + 
   geom_point(alpha=.5,size=3) +
   geom_errorbar(aes(ymin=z.S_cov - (1.96*sqrt(z.var.S_cov)),ymax=z.S_cov + (1.96*sqrt(z.var.S_cov)))) +
   geom_errorbarh(aes(xmin=z.N - (1.96*sqrt(z.var.N)),xmax=z.N + (1.96*sqrt(z.var.N)))) +
   geom_abline(intercept=0,slope=1,lty="dotted")

ggplot(data=ES_df.complete, aes(x=z.N,y=z.PIE,col=taxa)) + 
   geom_point(alpha=.5,size=3) +
   geom_errorbar(aes(ymin=z.PIE - (1.96*sqrt(z.var.PIE)),ymax=z.PIE + (1.96*sqrt(z.var.PIE)))) +
   geom_errorbarh(aes(xmin=z.N - (1.96*sqrt(z.var.N)),xmax=z.N + (1.96*sqrt(z.var.N)))) +
   geom_abline(intercept=0,slope=1,lty="dotted")

ggplot(data=ES_df.complete, aes(x=z.S,y=z.S_cov,col=taxa)) + 
   geom_point(alpha=.5,size=3) +
   geom_errorbar(aes(ymin=z.S_cov - (1.96*sqrt(z.var.S_cov)),ymax=z.S_cov + (1.96*sqrt(z.var.S_cov)))) +
   geom_errorbarh(aes(xmin=z.S - (1.96*sqrt(z.var.S)),xmax=z.S + (1.96*sqrt(z.var.S)))) +
   geom_abline(intercept=0,slope=1,lty="dotted")

ggplot(data=ES_df.complete, aes(x=z.S,y=z.PIE,col=taxa)) + 
   geom_point(alpha=.5,size=3) +
   geom_errorbar(aes(ymin=z.PIE - (1.96*sqrt(z.var.PIE)),ymax=z.PIE + (1.96*sqrt(z.var.PIE)))) +
   geom_errorbarh(aes(xmin=z.S - (1.96*sqrt(z.var.S)),xmax=z.S + (1.96*sqrt(z.var.S)))) +
   geom_abline(intercept=0,slope=1,lty="dotted")

ggplot(data=ES_df.complete, aes(x=z.S_cov,y=z.PIE,col=taxa)) + 
   geom_point(alpha=.5,size=3) +
   geom_errorbar(aes(ymin=z.PIE - (1.96*sqrt(z.var.PIE)),ymax=z.PIE + (1.96*sqrt(z.var.PIE)))) +
   geom_errorbarh(aes(xmin=z.S_cov - (1.96*sqrt(z.var.S_cov)),xmax=z.S + (1.96*sqrt(z.var.S_cov)))) +
   geom_abline(intercept=0,slope=1,lty="dotted")
dev.off()

############################################################################
### 1. Plot effect sizes in a 3-dim space
############################################################################
pdf(file=path2temp %+% "3D_Scatterplot_of_z.pdf")
with(ES_df.complete, scatter3D(x = z.N, y = z.S_cov, z = z.PIE, 
                               col = rainbow(length(levels(taxa))), 
                               pch = 16, cex = 1.5, alpha=0.6,
                               xlab = "Abundance N", ylab = "Species richness S", zlab = "ENS PIE", phi=5, ltheta=2, lphi=2,
                               #                       clab = c("Taxa"),
                               main = "", ticktype = "detailed",
                               colkey = F))
legend("bottomleft", paste(levels(ES_df.complete$taxa)), pch = 16, col = rainbow(length(levels(ES_df.complete$taxa))), cex=1, inset=c(-0.01,-0.03),bty="n")

with(ES_df.complete, scatter3D(x = z.N, z = z.S_cov, y = z.PIE, 
                               col = rainbow(length(levels(taxa))), 
                               pch = 16, cex = 1.5, alpha=0.6,
                               xlab = "Abundance N", zlab = "Species richness S", ylab = "ENS PIE", phi=5, ltheta=2, lphi=2,
                               #                       clab = c("Taxa"),
                               main = "", ticktype = "detailed",
                               colkey = F))
legend("bottomleft", paste(levels(ES_df.complete$taxa)), pch = 16, col = rainbow(length(levels(ES_df.complete$taxa))), cex=1, inset=c(-0.01,-0.03),bty="n")

with(ES_df.complete, scatter3D(z = z.N, y = z.S_cov, x = z.PIE, 
                               col = rainbow(length(levels(taxa))), 
                               pch = 16, cex = 1.5, alpha=0.6,
                               zlab = "Abundance N", ylab = "Species richness S", xlab = "ENS PIE", phi=5, ltheta=2, lphi=2,
                               #                       clab = c("Taxa"),
                               main = "", ticktype = "detailed",
                               colkey = F))
legend("bottomleft", paste(levels(ES_df.complete$taxa)), pch = 16, col = rainbow(length(levels(ES_df.complete$taxa))), cex=1, inset=c(-0.01,-0.03),bty="n")
dev.off()

############################################################################
### 2. forest plots
############################################################################
pdf(file=path2temp %+% "ForestPlots.pdf")
sapply(BDmetrics, function(j) forest(x=ES_df.complete[,"z." %+% j],vi=ES_df.complete[,"z.var." %+% j],slab=ES_df.complete$Study.ID,psize=1,main=paste(j),cex=.6))
dev.off()       

############################################################################
### 3. pairwise correlation of effect sizes based on complete observations
############################################################################

## put (absolute) correlations on the upper panels, with background color corresponing to the degree of correlation
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor=1, ...)
{
   usr <- par("usr"); on.exit(par(usr))
   par(usr = c(0, 1, 0, 1))
   r <- abs(cor(x, y, use="pairwise.complete.obs"))
   palette.vec <- colorRampPalette(colors=c("white","red"))(10)
   col.vec <- palette.vec[r*10]
   txt <- format(c(r, 0.123456789), digits = digits)[1]
   txt <- paste0(prefix, txt)
   if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
   rect(0,0,1,1,col=col.vec)
   text(0.5, 0.5, txt, cex = cex.cor)
}

pairs(ES_df.complete[,"z." %+% BDmetrics], labels=BDmetrics,lower.panel = panel.smooth, upper.panel = panel.cor, cex.labels=0.6)

png(file=path2temp %+% "PairwiseCorPlot.png", width=20,height=20,units="cm",res=400,type = "cairo-png")
pairs(ES_df.complete[,"z." %+% BDmetrics], labels=BDmetrics,lower.panel = panel.smooth, upper.panel = panel.cor, cex.labels=0.6)
dev.off()
