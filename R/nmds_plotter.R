nmds_plotter <- function(nmds,main=""){
  source("R/multiplot.R")
p <- 
  ggplot2::ggplot(data = nmds,aes(MDS1,MDS2,color=Sex,shape=Brood_Status)) +
  scale_colour_manual(values = c("hotpink","cornflowerblue"))   +
  geom_point(size=5) + geom_text(aes(label=Nest_Brood),size=2,colour="black") +
  coord_fixed()+
  labs(title =main, x = "", y = "")+
  theme_bw()+theme(
    plot.title=element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank())
p2 <- p + facet_grid(.~ Species) + theme(plot.title=element_blank(),legend.position = "none") 
multiplot(p,p2)
}