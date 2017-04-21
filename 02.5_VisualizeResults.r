load(path2temp %+% "Data4Analysis.Rdata") 
load(path2temp %+% "02.3_DataAnalysis_out.Rdata") # model_frag, model_frag_group, model_gradient

theme_habfrag <- function (base_size = 12, base_family = "", rel.text.size=1.5, legend.position = "right") {
   theme_grey(base_size = base_size, base_family = base_family) %+replace% 
      theme(axis.text = element_text(size = rel(rel.text.size-0.2)), 
            axis.ticks = element_line(colour = "black"), 
            axis.title = element_text(size = rel(rel.text.size)), 
            legend.text=element_text(size = rel(rel.text.size)),
            legend.title=element_text(size = rel(rel.text.size)),
            legend.key = element_rect(colour = "grey80"), 
            legend.position = legend.position,
            panel.background = element_rect(fill = "white", colour = NA), 
            panel.border = element_rect(fill = NA, colour = "grey50"), 
            panel.grid.major = element_line(colour = "grey90", size = 0.2), 
            panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
            strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2),
            strip.text = element_text(size=rel(rel.text.size)))
}

plot.func <- function(model,ylab,title){
   pred.se.list <- lapply(model,function(x) data.frame(b=x$b,se=x$se))
#   str(pred.se.list)
   pred.se.df <- bind_rows(pred.se.list,.id="model")
   pred.se.df$variable <- c("S", "D0_hat", "N_std", "ENS_pie", "D0_hat", "ENS_pie")
   ylim <- c(floor(10*min(pred.se.df$b-1.96*pred.se.df$se))/10,ceiling(10*max(pred.se.df$b+1.96*pred.se.df$se))/10)
   
   plot1 <- ggplot(pred.se.df[!(pred.se.df$model == "Common vs. rare"),]) +
      geom_point(aes(x=variable, y=b),size=4) +
      geom_errorbar(aes(x=variable, ymin=b-1.96*se,ymax=b+1.96*se),width=0.2,size=1.2) +
      geom_hline(yintercept=0,linetype="twodash", size=0.6) +
      scale_x_discrete("",limits=c("S", "D0_hat", "N_std", "ENS_pie")) +
      xlab("") + ylab(ylab) + ylim(ylim) +
      ggtitle(title) +
      theme_habfrag(rel.text.size=1.5)   
   
   plot2 <- ggplot(pred.se.df[(pred.se.df$model == "Common vs. rare"),]) +
      geom_point(aes(x=variable, y=b),size=4) +
      geom_errorbar(aes(x=variable, ymin=b-1.96*se,ymax=b+1.96*se),width=0.2,size=1.2) +
      geom_hline(yintercept=0,linetype="twodash", size=0.6) +
      scale_x_discrete("",limits=c("D0_hat", "ENS_pie")) +
      xlab("") + ylab("") + ylim(ylim) +
      ggtitle("") +      
      theme_habfrag(rel.text.size=1.5) +
      theme(axis.text.y=element_blank())  # remove y-axis labelling
      
   plot <- grid.arrange(plot1,plot2, nrow=1,ncol=2, widths=c(1,0.5))
   
   print(plot)
}

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func(model=model_frag, ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_group.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func(model_frag_group, ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_frag_gradient.png", width=20,height=10,units="cm",res=200,type = "cairo-png")
plot.func(model_gradient, ylab="Fishers' z",title="Gradient of fragmentation")
dev.off()

#####################################################################################
plot.func <- function(model,ylab,title){
   pred.se.list <- lapply(model,function(x) data.frame(b=x$b,se=x$se))
   #   str(pred.se.list)
   pred.se.df <- bind_rows(pred.se.list,.id="model")
   pred.se.df$variable <- rep(c("S", "D0_hat", "N_std", "ENS_pie", "D0_hat", "ENS_pie"),each=5)
   pred.se.df$variable[pred.se.df$model=="Common vs. rare"] <- rep(c("D0_hat", "ENS_pie"))
   pred.se.df$taxa <- c("invertebrates", "amphibians & reptiles", "birds", "mammals","plants") 
   pred.se.df$taxa[pred.se.df$model=="Common vs. rare"] <- rep(c("invertebrates", "amphibians & reptiles", "birds", "mammals","plants"), each=2)
   ylim <- c(floor(10*min(pred.se.df$b-1.96*pred.se.df$se))/10,ceiling(10*max(pred.se.df$b+1.96*pred.se.df$se))/10)
   
   pd <- position_dodge(width=0.6)
   
   plot1 <- ggplot(pred.se.df[!(pred.se.df$model == "Common vs. rare"),], aes(color=taxa)) +
      geom_point(position=pd,aes(x=variable, y=b),size=4) +
      geom_errorbar(position=pd,aes(x=variable, ymin=b-1.96*se,ymax=b+1.96*se),width=0.2,size=1.2) +
      geom_hline(yintercept=0,linetype="twodash", size=0.6) +
      scale_x_discrete("",limits=c("S", "D0_hat", "N_std", "ENS_pie")) +
      xlab("") + ylab(ylab) + ylim(ylim) +
      ggtitle(title) +
      guides(color=F) +
      theme_habfrag(rel.text.size=1.5)   
   
   plot2 <- ggplot(pred.se.df[(pred.se.df$model == "Common vs. rare"),], aes(color=taxa)) +
      geom_point(position=pd,aes(x=variable, y=b),size=4) +
      geom_errorbar(position=pd,aes(x=variable, ymin=b-1.96*se,ymax=b+1.96*se),width=0.2,size=1.2) +
      geom_hline(yintercept=0,linetype="twodash", size=0.6) +
      scale_x_discrete("",limits=c("D0_hat", "ENS_pie")) +
      xlab("") + ylab("") + ylim(ylim) +
      ggtitle("") +      
      theme_habfrag(rel.text.size=1.5) +
      theme(axis.text.y=element_blank())  # remove y-axis labelling
   
   plot <- grid.arrange(plot1,plot2, ncol=2, widths=c(1,1))
   
   print(plot)
}

png(file=path2temp %+% "ResultsPlots/ResultPlot_Taxafrag.png", width=30,height=10,units="cm",res=200,type = "cairo-png")
plot.func(modelTaxa_frag, ylab="Log(Response Ratio)",title="Smallest vs. largest fragment")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_Taxafrag_group.png", width=30,height=10,units="cm",res=200,type = "cairo-png")
plot.func(modelTaxa_frag_group, ylab="Log(Response Ratio)",title="Smallest vs. largest fragment group")
dev.off()

png(file=path2temp %+% "ResultsPlots/ResultPlot_Taxafrag_gradient.png", width=30,height=10,units="cm",res=200,type = "cairo-png")
plot.func(modelTaxa_gradient, ylab="Fishers' z",title="Gradient of fragmentation")
dev.off()

