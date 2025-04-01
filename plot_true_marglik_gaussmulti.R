rm(list=ls())
if(strsplit(getwd(),"/")[[1]][length(strsplit(getwd(),"/")[[1]])]!="thames_mixtures"){
  setwd("thames_mixtures")
}

library(ggplot2)


#df = read.csv('data/res_truemarglikgaussmulti.csv')
df = read.csv('data/res_truemarglikgaussmultiG15R5.csv') 
df_thames = as.data.frame(cbind(df[(df$thamesalg=="simple")&(df$estimate=="upper"),]$T,
                             df[(df$thamesalg=="simple")&(df$estimate=="estim"),]$value,
                             df[(df$thamesalg=="simple")&(df$estimate=="upper"),]$value,
                             df[(df$thamesalg=="simple")&(df$estimate=="lower"),]$value))
names(df_thames) = c("T","value","upper","lower")
df_bridge = as.data.frame(cbind(df[(df$thamesalg=="bridge")&(df$estimate=="upper"),]$T,
                                df[(df$thamesalg=="bridge")&(df$estimate=="estim"),]$value,
                                df[(df$thamesalg=="bridge")&(df$estimate=="upper"),]$value,
                                df[(df$thamesalg=="bridge")&(df$estimate=="lower"),]$value))
names(df_bridge) = c("T","value","upper","lower")
theme()

size_ggplot = 1.5
source("functions/theme.R")

annotations_G15R5 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("B"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1)) 

plot_G15R5 = ggplot(df[(df$thamesalg=="simple")&(df$estimate=="estim"),], aes(x=T,y=value)) +
  geom_point(color="blue",size=size_ggplot) +
  geom_point(color="red",data=df[(df$thamesalg=="bridge")&(df$estimate=="estim"),], aes(x=T,y=value),size = size_ggplot) + 
  labs(y=" ",title="G=15 components, dimension of d=5") +
  geom_abline(intercept = 0,slope=0) + labs(color='estimator') +theme(axis.text=element_text(size=12),
                                                                      axis.title=element_text(size=14,face="bold"),
                                                                      legend.text = element_text(size=12),
                                                                      legend.title = element_text(size=14,face="bold"),
                                                                      plot.title=element_text(size=13)) +
  scale_x_continuous(breaks = c(0, 50000, 100000, 150000, 200000), 
                     labels = c("0", "50000", "100000", "150000", "200000   ")) +
  geom_text(data=annotations_G15R5,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),size = 4)
#ggsave(filename="atelier/res_gaussmulti_G15R5.pdf",plot_G15R5,height=5,width=8)


df = read.csv('data/res_truemarglikgaussmulti.csv')
#df = read.csv('data/res_truemarglikgaussmultiG15R5.csv') 
df_thames = as.data.frame(cbind(df[(df$thamesalg=="simple")&(df$estimate=="upper"),]$T,
                                df[(df$thamesalg=="simple")&(df$estimate=="estim"),]$value,
                                df[(df$thamesalg=="simple")&(df$estimate=="upper"),]$value,
                                df[(df$thamesalg=="simple")&(df$estimate=="lower"),]$value))
names(df_thames) = c("T","value","upper","lower")
df_bridge = as.data.frame(cbind(df[(df$thamesalg=="bridge")&(df$estimate=="upper"),]$T,
                                df[(df$thamesalg=="bridge")&(df$estimate=="estim"),]$value,
                                df[(df$thamesalg=="bridge")&(df$estimate=="upper"),]$value,
                                df[(df$thamesalg=="bridge")&(df$estimate=="lower"),]$value))
names(df_bridge) = c("T","value","upper","lower")
theme()
source("functions/theme.R")

annotations_G5R6 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("A"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1)) 


plot_G5R6 = ggplot(df[(df$thamesalg=="simple")&(df$estimate=="estim"),], aes(x=T,y=value)) +
  geom_point(color="blue",size=size_ggplot) +
  geom_point(color="red",data=df[(df$thamesalg=="bridge")&(df$estimate=="estim"),], aes(x=T,y=value),size = size_ggplot) + 
  labs(y="error",title="G=5 components, dimension of d=6") +
  geom_abline(intercept = 0,slope=0) + labs(color='estimator') +theme(axis.text=element_text(size=12),
                                                                      axis.title=element_text(size=14,face="bold"),
                                                                      legend.text = element_text(size=12),
                                                                      legend.title = element_text(size=14,face="bold"),
                                                                      plot.title=element_text(size=13)) +
  scale_x_continuous(breaks = c(0, 50000, 100000, 150000, 200000), 
                     labels = c("0", "50000", "100000", "150000", "200000   ")) +
  geom_text(data=annotations_G5R6,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),size = 4)
#ggsave(filename="atelier/res_gaussmulti_G5R6.pdf",plot_G5R6,height=5,width=8)

library(gridExtra)
plot_gaussmulti_final = grid.arrange(plot_G5R6,plot_G15R5,ncol=2)
ggsave(filename="atelier/res_gaussmulti.pdf",plot_gaussmulti_final,height=2.5,width=8)

# Create a separate plot to act as the legend
data <- data.frame(
  x = c(1, 2),
  y = c(1, 2),
  group = c("A", "B")
)


legend_plot <- ggplot(data, aes(x = x, y = y, color = group)) +
  geom_point(data=data[1,],aes(color="THAMES"),size = size_ggplot) + 
  geom_point(data=data[2,],aes(color="bridge"),size = size_ggplot) +
  scale_color_manual(values = c("THAMES" = "blue","bridge" = "red")) +
  theme(
    legend.position = "top", # Legend positioning
    legend.title = element_text(size = 28),
    legend.text = element_text(size = 26),
  ) + labs(color = "estimators")


library(cowplot)

legend=get_plot_component(legend_plot,"guide-box-top")

ggsave(filename="atelier/res_gaussmulti_legend.pdf",legend,height=5,width=8)