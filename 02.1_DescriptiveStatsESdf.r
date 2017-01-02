div_df <- read.csv(path2temp %+% "DiversityData.csv", sep=";")
ES_df <- read.csv(file=path2temp %+% "ES_df.csv")

### number of studies with more than 3 fragments
ES_df <- subset(ES_df,n.fragment>3)
ES_df <- ES_df[complete.cases(ES_df[,"z." %+% names(div_df)[-(1:3)]]),]

### forest plots
pdf(path2temp %+% "ForestPlots.pdf")
sapply(names(div_df)[-(1:3)], function(j) forest(x=ES_df[,"z." %+% j],vi=ES_df[,"z.var." %+% j],slab=ES_df$Study.ID,psize=1,main=paste(j),cex=.6))
dev.off()       

### pairwise correlation of effect sizes based on complete observations

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

pairs(ES_df[,"z." %+% names(div_df)[-(1:3)]], labels=names(div_df)[-(1:3)],lower.panel = panel.smooth, upper.panel = panel.cor, cex.labels=0.6)

png(file=path2temp %+% "PairwiseCorPlot.png", width=20,height=20,units="cm",res=400,type = "cairo-png")
pairs(ES_df[,"z." %+% names(div_df)[-(1:3)]], labels=names(div_df)[-(1:3)],lower.panel = panel.smooth, upper.panel = panel.cor, cex.labels=0.6)
dev.off()
