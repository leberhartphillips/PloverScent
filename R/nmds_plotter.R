nmds_plotter <- function(nmds,main=""){
  source("R/multiplot.R")
p <- 
  ggplot2::ggplot(data = nmds,aes(MDS1,MDS2,color=as.factor(Sex),shape=Species)) +
  geom_point(size=7) + geom_text(aes(label=GC_Sample),size=2,colour="black") +
  scale_x_continuous(limits = c(-1.25,1.25)) +
  scale_y_continuous(limits = c(-1.25,1.25)) +
  scale_shape_manual(values=c(15,16,17)) 
  labs(title =main, x = "MDS1", y = "MDS2")+
  theme_bw()+theme(
    plot.title=element_text(face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.text.x  = element_blank(), 
    axis.title.y = element_text(size = 16),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank())
# geom_rect(colour="black",fill="khaki1",aes(xmin=1,xmax=1.2,ymin=0.95,ymax=1.25),alpha=0.2) +
# annotate("text",x=1.1,y=1.1,label=paste(c("Adonis",R2,Pval),collapse = "\n"),col="Black",size=5)
p2 <- p + facet_grid(.~ Species) + theme(plot.title=element_blank(),legend.position = "none")
multiplot(p,p2)
}