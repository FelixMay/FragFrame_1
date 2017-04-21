load(path2temp %+% "Data4Analysis.Rdata") 
ls()

BDmetrics <- c("S","D0_hat","N_std","ENS_pie")
############################################################################
### 1. Forest plots
############################################################################

############################################################################
### for each BD metric separately
pdf(file=path2temp %+% "ForestPlots_z.pdf")
sapply(BDmetrics, function(j) forest(x=ES_df.complete[,"ES." %+% j],vi=ES_df.complete[,"ES.var." %+% j],slab=ES_df.complete$Study.ID,psize=1,main=paste(j),cex=.6))
dev.off()       

############################################################################
### combined forest plot for complete dataset
png(path2temp %+% "ForestPlots_ES_df.complete_extremes.png",height=40,width=25,units="cm",res=400)
g_legend<-function(a.gplot){ 
   tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
   legend <- tmp$grobs[[leg]] 
   return(legend)} 

### save legend
plot1 <- ggplot(data=ES_df.complete_long,aes(y=Case.ID,yend=Case.ID,x=value-1.96*sqrt(ES.var),xend=value+1.96*sqrt(ES.var),color=ES)) + 
   geom_segment() +
   geom_point(aes(x=value)) +
   scale_color_brewer("", palette="Set1",breaks=c("ES.S","ES.D0_hat","ES.N_std","ES.ENS_pie"), labels=c("S","D0_hat","N_std","ENS_PIE"))

# # if legend from plot1 is drawn
legend <- g_legend(plot1)
# png(path2temp %+% "legend.png",height=10,width=10,units="cm",res=200)
# grid.draw(legend)
# dev.off()

pd <- position_dodge(width=0.6)
plot2 <- ggplot(data=ES_df.complete_long, aes(x=Case.ID,y=value,ymin=value-1.96*sqrt(ES.var),ymax=value+1.96*sqrt(ES.var),color=ES)) + 
   geom_pointrange(position=pd, cex=0.4,show.legend = F) + 
   geom_hline(aes(yintercept=0), linetype="twodash",size=0.6) +
   xlab("") + ylab("") +
   coord_flip() +
   scale_color_brewer("", palette="Set1",breaks=c("ES.S","ES.D0_hat","ES.N_std","ES.ENS_pie"), labels=c("S","D0_hat","N_std","ENS_PIE"))

grid.arrange(plot2,legend, ncol=2, widths=c(15,4))
dev.off()

############################################################################
### combined forest plot for complete dataset without extremes of Goodman and Zartman
png(path2temp %+% "ForestPlots_ES_df.complete_wo_extremes.png",height=40,width=25,units="cm",res=400)
g_legend<-function(a.gplot){ 
   tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
   legend <- tmp$grobs[[leg]] 
   return(legend)} 

### save legend
plot1 <- ggplot(data=ES_df.complete_long,aes(y=Case.ID,yend=Case.ID,x=value-1.96*sqrt(ES.var),xend=value+1.96*sqrt(ES.var),color=ES)) + 
         geom_segment() +
         geom_point(aes(x=value)) +
         scale_color_brewer("", palette="Set1",breaks=c("ES.S","ES.D0_hat","ES.N_std","ES.ENS_pie"), labels=c("S","D0_hat","N_std","ENS_PIE"))

# # if legend from plot1 is drawn
legend <- g_legend(plot1)
# png(path2temp %+% "legend.png",height=10,width=10,units="cm",res=200)
# grid.draw(legend)
# dev.off()

pd <- position_dodge(width=0.6)
plot2 <- ggplot(data=ES_df.complete_long, aes(x=Case.ID,y=value,ymin=value-1.96*sqrt(ES.var),ymax=value+1.96*sqrt(ES.var),color=ES)) + 
   geom_pointrange(position=pd, cex=0.4,show.legend = F) + 
   geom_hline(aes(yintercept=0), linetype="twodash",size=0.6) +
   ylim(-5,5) +
   xlab("") + ylab("") +
   coord_flip() +
   scale_color_brewer("", palette="Set1",breaks=c("ES.S","ES.D0_hat","ES.N_std","ES.ENS_pie"), labels=c("S","D0_hat","N_std","ENS_PIE"))

grid.arrange(plot2,legend, ncol=2, widths=c(15,4))
dev.off()

############################################################################
### 2. pairwise correlation of effect sizes based on complete observations
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
panel.smooth <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                          cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
   points(x, y, pch = pch, col = col, bg = bg, cex = cex)
   ok <- is.finite(x) & is.finite(y)
   if (any(ok)) 
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)
   abline(0,1,lty="dotted")
}

png(file=path2temp %+% "PairwiseCorPlots/PairwiseCorPlot_frag.png", width=20,height=20,units="cm",res=400,type = "cairo-png")
pairs(ES_frag_df.complete[,"ES." %+% BDmetrics], labels=BDmetrics,lower.panel = panel.smooth, upper.panel = panel.cor)
dev.off()
png(file=path2temp %+% "PairwiseCorPlots/PairwiseCorPlot_frag_group.png", width=20,height=20,units="cm",res=400,type = "cairo-png")
pairs(ES_frag_group_df.complete[,"ES." %+% BDmetrics], labels=BDmetrics,lower.panel = panel.smooth, upper.panel = panel.cor)
dev.off()
png(file=path2temp %+% "PairwiseCorPlots/PairwiseCorPlot_gradient.png", width=20,height=20,units="cm",res=400,type = "cairo-png")
pairs(ES_df.complete[,"ES." %+% BDmetrics], labels=BDmetrics,lower.panel = panel.smooth, upper.panel = panel.cor)
dev.off()

############################################################################
### 3. Histograms of predictor variables
############################################################################

for(col in c("taxa","country", "continent", "biome", "fragment.biome","matrix.biome", "fragment.veg", "matrix.veg")){
   p <- ggplot(data=meta_df) + 
      geom_histogram(aes(x=meta_df[,col]), size=0.4, binwidth=0.2,stat="count") + 
      labs(x="",y="") +
      ggtitle(paste(col)) + 
      theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)),plot.title=element_text(size = rel(2)) , axis.text.x=element_text(angle=45,vjust = 1, hjust=1),legend.text=element_text(size = rel(2)),legend.title=element_text(size = rel(2)))
   #,axis.ticks.length=unit(.4,"cm")
   print(p)
   ggsave(p, file = path2temp %+% "Histograms/Histogram_meta_df_" %+% col %+% ".png", width = 20, height = 8, type = "cairo-png")
   
}


#### RESTERAMPE
############################################################################
### 1. Plot effect sizes in a 2-dim space
############################################################################
# pdf(file=path2temp %+% "2D_Scatterplot_of_z.pdf")
# ggplot(data=ES_df.complete, aes(x=z.N,y=z.S,col=taxa)) + 
#    geom_point(alpha=.5,size=3) +
#    geom_errorbar(aes(ymin=z.S - (1.96*sqrt(z.var.S)),ymax=z.S + (1.96*sqrt(z.var.S)))) +
#    geom_errorbarh(aes(xmin=z.N - (1.96*sqrt(z.var.N)),xmax=z.N + (1.96*sqrt(z.var.N)))) +
#    geom_abline(intercept=0,slope=1,lty="dotted")
# 
# ggplot(data=ES_df.complete, aes(x=z.N,y=z.ENS_pie,col=taxa)) + 
#    geom_point(alpha=.5,size=3) +
#    geom_errorbar(aes(ymin=z.ENS_pie - (1.96*sqrt(z.var.ENS_pie)),ymax=z.ENS_pie + (1.96*sqrt(z.var.ENS_pie)))) +
#    geom_errorbarh(aes(xmin=z.N - (1.96*sqrt(z.var.N)),xmax=z.N + (1.96*sqrt(z.var.N)))) +
#    geom_abline(intercept=0,slope=1,lty="dotted")
# 
# ggplot(data=ES_df.complete, aes(x=z.S,y=z.ENS_pie,col=taxa)) + 
#    geom_point(alpha=.5,size=3) +
#    geom_errorbar(aes(ymin=z.ENS_pie - (1.96*sqrt(z.var.ENS_pie)),ymax=z.ENS_pie + (1.96*sqrt(z.var.ENS_pie)))) +
#    geom_errorbarh(aes(xmin=z.S - (1.96*sqrt(z.var.S)),xmax=z.S + (1.96*sqrt(z.var.S)))) +
#    geom_abline(intercept=0,slope=1,lty="dotted")
# 
# dev.off()

# ############################################################################
# ### 1. Plot effect sizes in a 3-dim space
# ############################################################################
# pdf(file=path2temp %+% "3D_Scatterplot_of_z.pdf")
# with(ES_df.complete, scatter3D(x = z.N, y = z.S, z = z.ENS_pie, 
#                                col = rainbow(length(levels(taxa))), 
#                                pch = 16, cex = 1.5, alpha=0.6,
#                                xlab = "Abundance N", ylab = "Species richness S", zlab = "ENS PIE", phi=5, ltheta=2, lphi=2,
#                                #                       clab = c("Taxa"),
#                                main = "", ticktype = "detailed",
#                                colkey = F))
# legend("bottomleft", paste(levels(ES_df.complete$taxa)), pch = 16, col = rainbow(length(levels(ES_df.complete$taxa))), cex=1, inset=c(-0.01,-0.03),bty="n")
# 
# with(ES_df.complete, scatter3D(x = z.N, z = z.S, y = z.ENS_pie, 
#                                col = rainbow(length(levels(taxa))), 
#                                pch = 16, cex = 1.5, alpha=0.6,
#                                xlab = "Abundance N", zlab = "Species richness S", ylab = "ENS PIE", phi=5, ltheta=2, lphi=2,
#                                #                       clab = c("Taxa"),
#                                main = "", ticktype = "detailed",
#                                colkey = F))
# legend("bottomleft", paste(levels(ES_df.complete$taxa)), pch = 16, col = rainbow(length(levels(ES_df.complete$taxa))), cex=1, inset=c(-0.01,-0.03),bty="n")
# 
# with(ES_df.complete, scatter3D(z = z.N, y = z.S, x = z.ENS_pie, 
#                                col = rainbow(length(levels(taxa))), 
#                                pch = 16, cex = 1.5, alpha=0.6,
#                                zlab = "Abundance N", ylab = "Species richness S", xlab = "ENS PIE", phi=5, ltheta=2, lphi=2,
#                                #                       clab = c("Taxa"),
#                                main = "", ticktype = "detailed",
#                                colkey = F))
# legend("bottomleft", paste(levels(ES_df.complete$taxa)), pch = 16, col = rainbow(length(levels(ES_df.complete$taxa))), cex=1, inset=c(-0.01,-0.03),bty="n")
# dev.off()

