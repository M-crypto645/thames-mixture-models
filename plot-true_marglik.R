# simulate from the Gaussian mixture model with known proportions and variances
# and compute the marginal likelihood estimators
rm(list=ls())
if(strsplit(getwd(),"/")[[1]][length(strsplit(getwd(),"/")[[1]])]!="thames_mixtures"){
  setwd("thames_mixtures")
}
library(ggplot2)
library(reshape2)
source("functions/theme.R")

df_G2 = read.csv('data/res_sim1_G2.csv')
df_overfitting = read.csv('data/res_sim1_overfitting.csv')
df_underfitting = read.csv('data/res_sim1_underfitting.csv')
df_G3 = read.csv('data/res_sim1_G3.csv')

get_plot = function(df){
  # reorder for ggplot
  df$estimator = df$thamesalg
  df$estimator <- factor(df$estimator , levels=c("mc", "bridge", "standard", "simple"))
  
  ggplot(df, aes(y=value,x=estimator)) + 
    geom_boxplot() + 
    geom_hline(yintercept=0,color="red") +
    theme(text=element_text(size=48))+
    scale_x_discrete(labels = c("mc","bridge",
                                "Reichl","THAMES")) + 
    labs(y="error")
}

get_plot_nomc = function(df){
  df$rho = 0
  df$rho = sapply(df$sampler,function(s) grepl("rho1",s))*1+sapply(df$sampler,function(s) grepl("rho.5",s))*0.5
  df$rho = c("0",df$rho)[-1]
  df$estimator = df$thamesalg
  df$estimator <- factor(df$estimator , levels=c("mc", "bridge", "standard", "simple"))
  
  ggplot(df, aes(y=value,x=rho)) + 
    geom_boxplot(aes(y=value,x=rho,fill=estimator),outlier.shape = NA) + 
    geom_hline(yintercept=0,color="red") +
    labs(y="error",x=expression(rho),main="true model") +
    geom_vline(xintercept = 1.5) + geom_vline(xintercept = 2.5) +
    scale_fill_manual(values = c("brown","grey","beige"),
                      labels =  expression("bridge","Reichl","THAMES"), 
                      name="estimators") + theme(legend.position = "none")
}

# reorder for ggplot
# df = df_G2[(df_G2$ellipsealg =="standard")&(df_G2$estimate=="estim")&(df_G3$thamesalg!="mc"),]
# df = df_G3[(df_G3$ellipsealg =="standard")&(df_G3$estimate=="estim")&(df_G3$thamesalg!="mc"),]
# df = df_overfitting[(df_overfitting$ellipsealg =="standard")&(df_overfitting$estimate=="estim")&(df_overfitting$thamesalg!="mc"),]

annotations_G2G2 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("A"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))
annotations_G2G3 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("B"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))
annotations_G3G2 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("C"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))

library(gridExtra)
plot_final = grid.arrange(get_plot_nomc(df_G2[(df_G2$ellipsealg =="standard")&(df_G2$estimate=="estim")&(df_G2$thamesalg!="mc"),]) + 
               labs(title="true model\n(fitted G=2, true G=2)") + 
                 geom_text(data=annotations_G2G2,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),size = 5),
             get_plot_nomc(df_underfitting[(df_underfitting$ellipsealg =="standard")&(df_underfitting$estimate=="estim")&(df_underfitting$thamesalg!="mc"),]) + 
               labs(title="underfitting\n(fitted G=2, true G=3)",y=" ") + 
               geom_text(data=annotations_G2G3,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),size = 5),
             get_plot_nomc(df_overfitting[(df_overfitting$ellipsealg =="standard")&(df_overfitting$estimate=="estim")&(df_overfitting$thamesalg!="mc"),]) + 
               labs(title="overfitting\n(fitted G=3, true G=2)",y=" ") +
               geom_text(data=annotations_G3G2,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),size = 5),ncol=3)
ggsave("atelier/res_gaussuniv.pdf",plot_final,height=2.5,width=8)

library(cowplot)

legend=get_plot_component(get_plot_nomc(df_G2[(df_G2$ellipsealg =="standard")&(df_G2$estimate=="estim")&(df_G3$thamesalg!="mc"),])+theme(legend.position = "top"),"guide-box-top")
ggsave(filename="atelier/res_gaussuniv_legend.pdf",legend,height=5,width=8)


### OLD PLOTS ###

# # true model (G=2)
# ggsave(get_plot(df_G2[(df_G2$ellipsealg =="standard")&(df_G2$estimate=="estim")&(df_G2$sampler=="truemodel_G2_rho0"),]),
#        width=10.6,height=8.14,file="atelier/plot_g2_true_model_rho_0.pdf")
# ggsave(get_plot(df_G2[(df_G2$ellipsealg =="standard")&(df_G2$estimate=="estim")&(df_G2$sampler=="truemodel_G2_rho.5"),]),
#        width=10.6,height=8.14,file="atelier/plot_g2_true_model_rho_.5.pdf")
# ggsave(get_plot(df_G2[(df_G2$ellipsealg =="standard")&(df_G2$estimate=="estim")&(df_G2$sampler=="truemodel_G2_rho1"),]),
#        width=10.6,height=8.14,file="atelier/plot_g2_true_model_rho_1.pdf")
# 
# # overfitting
# ggsave(get_plot(df_overfitting[(df_overfitting$ellipsealg =="standard")&(df_overfitting$estimate=="estim")&(df_overfitting$sampler=="overfitting_rho0"),]),
#        width=10.6,height=8.14,file="atelier/plot_overfitting_rho_0.pdf")
# ggsave(get_plot(df_overfitting[(df_overfitting$ellipsealg =="standard")&(df_overfitting$estimate=="estim")&(df_overfitting$sampler=="overfitting_rho.5"),]),
#        width=10.6,height=8.14,file="atelier/plot_overfitting_rho_.5.pdf")
# ggsave(get_plot(df_overfitting[(df_overfitting$ellipsealg =="standard")&(df_overfitting$estimate=="estim")&(df_overfitting$sampler=="overfitting_rho1"),]),
#        width=10.6,height=8.14,file="atelier/plot_overfitting_rho_1.pdf")
# 
# # underfitting
# ggsave(get_plot(df_underfitting[(df_underfitting$ellipsealg =="standard")&(df_underfitting$estimate=="estim")&(df_underfitting$sampler=="underfitting_rho0"),]),
#        width=10.6,height=8.14,file="atelier/plot_underfitting_rho_0.pdf")
# ggsave(get_plot(df_underfitting[(df_underfitting$ellipsealg =="standard")&(df_underfitting$estimate=="estim")&(df_underfitting$sampler=="underfitting_rho.5"),]),
#        width=10.6,height=8.14,file="atelier/plot_underfitting_rho_.5.pdf")
# ggsave(get_plot(df_underfitting[(df_underfitting$ellipsealg =="standard")&(df_underfitting$estimate=="estim")&(df_underfitting$sampler=="underfitting_rho1"),]),
#        width=10.6,height=8.14,file="atelier/plot_underfitting_rho_1.pdf")
# 
# # true model (G=3)
# ggsave(get_plot(df_G3[(df_G3$ellipsealg =="standard")&(df_G3$estimate=="estim")&(df_G3$sampler=="truemodel_G3_rho0"),]),
#        width=10.6,height=8.14,file="atelier/plot_g3_true_model_rho_0.pdf")
# ggsave(get_plot(df_G3[(df_G3$ellipsealg =="standard")&(df_G3$estimate=="estim")&(df_G3$sampler=="truemodel_G3_rho.5"),]),
#        width=10.6,height=8.14,file="atelier/plot_g3_true_model_rho_.5.pdf")
# ggsave(get_plot(df_G3[(df_G3$ellipsealg =="standard")&(df_G3$estimate=="estim")&(df_G3$sampler=="truemodel_G3_rho1"),]),
#        width=10.6,height=8.14,file="atelier/plot_g3_true_model_rho_1.pdf")