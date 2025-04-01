# Fit the THAMES, fitted to the galaxies, enzyme, and acidity datasets
rm(list=ls())

if(strsplit(getwd(),"/")[[1]][length(strsplit(getwd(),"/")[[1]])]!="thames_mixtures"){
  setwd("thames_mixtures")
}

library(ggplot2)
library(gridExtra)

# load(file=paste0('data/res_acidity','.Rda'))
# load(file=paste0('data/res_enzyme','.Rda'))
# load(file=paste0('data/res_galaxies','.Rda'))
df_acidity = read.csv('data/res_acidity.csv')
df_enzymes = read.csv('data/res_enzyme.csv')
df_galaxies = read.csv('data/res_galaxies.csv')

# fix 'continuous x aesthetic' bug
df_galaxies$G = c("",df_galaxies$G)[-1]
df_enzymes$G = c("",df_enzymes$G)[-1]
df_acidity$G = c("",df_acidity$G)[-1]

### BEGIN plot used in the main documents ###

# Celeux et al. (2018) galaxy data results
lml_thames_galaxies <- data.frame(G = 2:6, 
                                  lml = sapply(2:6, function(s) median(df_galaxies[(df_galaxies$sampler=="JAGS")&
                                                                  (df_galaxies$relabalg=="ECR")&
                                                                  (df_galaxies$ellipsealg=="standard")&
                                                                  (df_galaxies$thamesalg=="permutations")&
                                                                  (df_galaxies$estimate=="estim")&
                                                                    df_galaxies$G==paste0(s),]$value)),
                                  lower = sapply(2:6, function(s) quantile(df_galaxies[(df_galaxies$sampler=="JAGS")&
                                                                                       (df_galaxies$relabalg=="ECR")&
                                                                                       (df_galaxies$ellipsealg=="standard")&
                                                                                       (df_galaxies$thamesalg=="permutations")&
                                                                                       (df_galaxies$estimate=="estim")&
                                                                                       df_galaxies$G==paste0(s),]$value,0.025)),
                                  upper = sapply(2:6, function(s) quantile(df_galaxies[(df_galaxies$sampler=="JAGS")&
                                                                                         (df_galaxies$relabalg=="ECR")&
                                                                                         (df_galaxies$ellipsealg=="standard")&
                                                                                         (df_galaxies$thamesalg=="permutations")&
                                                                                         (df_galaxies$estimate=="estim")&
                                                                                         df_galaxies$G==paste0(s),]$value,1-0.025)))
  
lml_true_galaxies <- data.frame(G = 2:6, lml = c(-235.37,-226.85,-226.11,-225.83,-225.75),
                                y_higher=-222.091, y_lower=-240)
lml_ISF_galaxies  <- data.frame(G = 2:6, lml = c(-235.37,-226.85,-226.11,-225.93,-226.10),
                               y_higher=-222.091, y_lower=-240)
lml_ISR_galaxies  <- data.frame(G = 2:6, lml = c(-235.37,-226.85,-226.11,-225.93,-224.91),
                               y_higher=-222.091, y_lower=-240)
lml_BSR_galaxies  <- data.frame(G = 2:6, lml = c(-235.37,-226.85,-226.27,-227.20,-230.00),
                               y_higher=-222.091, y_lower=-240)
lml_RIF_galaxies  <- data.frame(G = 2:6, lml = c(-235.37,-227.62,-226.35,-226.18,-229.57),
                               y_higher=-222.091, y_lower=-240)
lml_RIR_galaxies  <- data.frame(G = 2:6, lml = c(-235.37,-227.03,-227.54,-230.00,-235.93),
                               y_higher=-222.091, y_lower=-240)

# Celeux et al. (2018) enzyme data results
#lml_thames_enzyme <- data.frame(G = 2:6, lml = )

lml_thames_enzyme <- data.frame(G = 2:6, lml = sapply(2:6, function(s) median(df_enzymes[(df_enzymes$sampler=="JAGS")&
                                                                                           (df_enzymes$relabalg=="ECR")&
                                                                                           (df_enzymes$ellipsealg=="standard")&
                                                                                           (df_enzymes$thamesalg=="permutations")&
                                                                                           (df_enzymes$estimate=="estim")&
                                                                                           (df_enzymes$G==paste0(s)),]$value)),
                                lower = sapply(2:6, function(s) quantile(df_enzymes[(df_enzymes$sampler=="JAGS")&
                                                                                    (df_enzymes$relabalg=="ECR")&
                                                                                    (df_enzymes$ellipsealg=="standard")&
                                                                                    (df_enzymes$thamesalg=="permutations")&
                                                                                    (df_enzymes$estimate=="estim")&
                                                                                    (df_enzymes$G==paste0(s)),]$value,0.025)),
                                upper = sapply(2:6, function(s) quantile(df_enzymes[(df_enzymes$sampler=="JAGS")&
                                                                                      (df_enzymes$relabalg=="ECR")&
                                                                                      (df_enzymes$ellipsealg=="standard")&
                                                                                      (df_enzymes$thamesalg=="permutations")&
                                                                                      (df_enzymes$estimate=="estim")&
                                                                                      (df_enzymes$G==paste0(s)),]$value,1-0.025)))

lml_true_enzyme <- data.frame(G = 2:6, lml = c(-76.5,-74.2,-74.3,-75,-76.7),
                              y_higher=-73, y_lower=-90)
lml_ISF_enzyme <- data.frame(G = 2:6, lml = c(-76.5,-74.2,-74.3,-75.16,-77.37),
                              y_higher=-73, y_lower=-90)
lml_ISR_enzyme <- data.frame(G = 2:6, lml = c(-76.5,-73.89,-73.89,-75.25,-79.66),
                             y_higher=-73, y_lower=-90)
lml_BSR_enzyme <- data.frame(G = 2:6, lml = c(-76.5,-74.2,-74.74,-79.32,-84.49),
                             y_higher=-73, y_lower=-90)
lml_RIF_enzyme <- data.frame(G = 2:6, lml = c(-76.5,-74.2,-75.08,-79.15,-82.28),
                             y_higher=-73, y_lower=-90)
lml_RIR_enzyme <- data.frame(G = 2:6, lml = c(-76.5,-74.2,-77.71,-83.72,-89.06),
                             y_higher=-73, y_lower=-90)

# Celeux et al. (2018) acidity data results
lml_thames_acidity <- data.frame(G=2:6, lml = sapply(2:6, function(s) median(df_acidity[(df_acidity$sampler=="JAGS")&
                                                                                          (df_acidity$relabalg=="ECR")&
                                                                                          (df_acidity$ellipsealg=="standard")&
                                                                                          (df_acidity$thamesalg=="permutations")&
                                                                                          (df_acidity$estimate=="estim")&(df_acidity$G==paste0(s)),]$value)))

lml_true_acidity <- data.frame(G = 2:6, lml = c(-199.4,-198.2,-198.3,-199.0,-200.1),
                               y_lower=-216.988,y_higher=-194.064)
lml_ISF_acidity <- data.frame(G = 2:6, lml = c(-199.4,-198.2,-198.26,-200.0,-200.54),
                               y_lower=-216.988,y_higher=-194.064)
lml_ISR_acidity <- data.frame(G = 2:6, lml = c(-199.4,-198.2,-198.58,-198.26,-197.17),
                               y_lower=-216.988,y_higher=-194.064)
lml_BSR_acidity <- data.frame(G = 2:6, lml = c(-199.4,-198.2,-199.02,-203.47,-207.5),
                               y_lower=-216.988,y_higher=-194.064)
lml_RIF_acidity <- data.frame(G = 2:6, lml = c(-199.4,-198.2,-199.45,-201.95,-205.65),
                               y_lower=-216.988,y_higher=-194.064)
lml_RIR_acidity <- data.frame(G = 2:6, lml = c(-199.4,-198.36,-201.95,-207.28,-213.69),
                               y_lower=-216.988,y_higher=-194.064)

# graphical parameters
size_ggplot = 2
# colors = c("green3", "red", "blue", "pink", "turquoise", "brown", "black")
#colors = c("green3", "green3", "red", "red", "red", "red", "black")
# colors = c("blue", "blue", "red", "red", "red", "red", "black")

colors = c("brown","grey","pink","pink3","darkorange2","orange2","black")
library(viridis)
library(RColorBrewer)
RColorBrewer::brewer.pal(5,"YlOrRd")
#colors = c("blue", "#1C9099",
#           "black", "red", "orange2", "brown", "cyan")
#colors = c("blue","cyan3",viridis::magma(12)[7:10],"black")
#c("#000004FF" "#2D1160FF" "#721F81FF" "#B63679FF" "#F1605DFF", "#FEAF77FF", "#FCFDBFFF")
#shapes = 0:6
shapes = c(rep(16,6),15)
startx = 1.6
mindist = 0.1325

source("functions/theme.R")
theme_set(my_theme)
plotgalaxies = ggplot(df_galaxies[(df_galaxies$sampler=="JAGS")&(df_galaxies$relabalg=="ECR")&(df_galaxies$ellipsealg=="standard")&(df_galaxies$thamesalg=="permutations")&(df_galaxies$estimate=="estim"),],
                      aes(x=G,y=value)) + annotate("point", x = startx, y = lml_true_galaxies$lml[1], colour = colors[1],shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 1 + startx, y = lml_true_galaxies$lml[2], colour = colors[1], shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 2 + startx, y = lml_true_galaxies$lml[3], colour = colors[1], shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 3 + startx, y = lml_true_galaxies$lml[4], colour = colors[1], shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 4 + startx, y = lml_true_galaxies$lml[5], colour = colors[1], shape=shapes[1], size=size_ggplot) +
  annotate("point", x = startx + mindist, y = lml_ISF_galaxies$lml[1], colour = colors[2], shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 1 + startx + mindist, y = lml_ISF_galaxies$lml[2], colour = colors[2], shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 2 + startx + mindist, y = lml_ISF_galaxies$lml[3], colour = colors[2], shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 3 + startx + mindist, y = lml_ISF_galaxies$lml[4], colour = colors[2], shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 4 + startx + mindist, y = lml_ISF_galaxies$lml[5], colour = colors[2], shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 0 + startx + 3*mindist, y = lml_ISR_galaxies$lml[1], colour = colors[3], shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 1 + startx + 3*mindist, y = lml_ISR_galaxies$lml[2], colour = colors[3], shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 2 + startx + 3*mindist, y = lml_ISR_galaxies$lml[3], colour = colors[3], shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 3 + startx + 3*mindist, y = lml_ISR_galaxies$lml[4], colour = colors[3], shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 4 + startx + 3*mindist, y = lml_ISR_galaxies$lml[5], colour = colors[3], shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 0 + startx + 4*mindist, y = lml_BSR_galaxies$lml[1], colour = colors[4], shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 1 + startx + 4*mindist, y = lml_BSR_galaxies$lml[2], colour = colors[4], shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 2 + startx + 4*mindist, y = lml_BSR_galaxies$lml[3], colour = colors[4], shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 3 + startx + 4*mindist, y = lml_BSR_galaxies$lml[4], colour = colors[4], shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 4 + startx + 4*mindist, y = lml_BSR_galaxies$lml[5], colour = colors[4], shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 0 + startx + 5*mindist, y = lml_RIF_galaxies$lml[1], colour = colors[5], shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 1 + startx + 5*mindist, y = lml_RIF_galaxies$lml[2], colour = colors[5], shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 2 + startx + 5*mindist, y = lml_RIF_galaxies$lml[3], colour = colors[5], shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 3 + startx + 5*mindist, y = lml_RIF_galaxies$lml[4], colour = colors[5], shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 4 + startx + 5*mindist, y = lml_RIF_galaxies$lml[5], colour = colors[5], shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 0 + startx + 6*mindist, y = lml_RIR_galaxies$lml[1], colour = colors[6], shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 1 + startx + 6*mindist, y = lml_RIR_galaxies$lml[2], colour = colors[6], shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 2 + startx + 6*mindist, y = lml_RIR_galaxies$lml[3], colour = colors[6], shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 3 + startx + 6*mindist, y = lml_RIR_galaxies$lml[4], colour = colors[6], shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 4 + startx + 6*mindist, y = lml_RIR_galaxies$lml[5], colour = colors[6], shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 0 + startx + 2*mindist, y = lml_thames_galaxies$lml[1], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  annotate("point", x = 1 + startx + 2*mindist, y = lml_thames_galaxies$lml[2], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  annotate("point", x = 2 + startx + 2*mindist, y = lml_thames_galaxies$lml[3], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  annotate("point", x = 3 + startx + 2*mindist, y = lml_thames_galaxies$lml[4], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  annotate("point", x = 4 + startx + 2*mindist, y = lml_thames_galaxies$lml[5], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  geom_vline(xintercept=2.5) +
  geom_vline(xintercept=3.5) +
  geom_vline(xintercept=4.5) +
  geom_vline(xintercept=5.5) +
  geom_vline(xintercept=6.5) + labs(title="galaxy dataset",y="log marginal likelihood") + coord_cartesian(ylim=c(-237,-222),xlim = c(1.7,6.3))# xlim(0.6,6.6)

annotations_galaxies <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("A"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1)) 

plotgalaxies = plotgalaxies + theme(axis.title=element_text(size=10)) + 
  geom_text(data=annotations_galaxies,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),size = 5)
plotenzyme = ggplot(df_enzymes[(df_enzymes$sampler=="JAGS")&(df_enzymes$relabalg=="ECR")&(df_enzymes$ellipsealg=="standard")&(df_enzymes$thamesalg=="permutations")&(df_enzymes$estimate=="estim"),],
                    aes(x=G,y=value)) + annotate("point", x = 0 + startx + 0*mindist, y = lml_true_enzyme$lml[1], colour = colors[1],shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 1 + startx + 0*mindist, y = lml_true_enzyme$lml[2], colour = colors[1],shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 2 + startx + 0*mindist, y = lml_true_enzyme$lml[3], colour = colors[1],shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 3 + startx + 0*mindist, y = lml_true_enzyme$lml[4], colour = colors[1],shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 4 + startx + 0*mindist, y = lml_true_enzyme$lml[5], colour = colors[1],shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 0 + startx + 1*mindist, y = lml_ISF_enzyme$lml[1], colour = colors[2],shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 1 + startx + 1*mindist, y = lml_ISF_enzyme$lml[2], colour = colors[2],shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 2 + startx + 1*mindist, y = lml_ISF_enzyme$lml[3], colour = colors[2],shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 3 + startx + 1*mindist, y = lml_ISF_enzyme$lml[4], colour = colors[2],shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 4 + startx + 1*mindist, y = lml_ISF_enzyme$lml[5], colour = colors[2],shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 0 + startx + 3*mindist, y = lml_ISR_enzyme$lml[1], colour = colors[3],shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 1 + startx + 3*mindist, y = lml_ISR_enzyme$lml[2], colour = colors[3],shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 2 + startx + 3*mindist, y = lml_ISR_enzyme$lml[3], colour = colors[3],shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 3 + startx + 3*mindist, y = lml_ISR_enzyme$lml[4], colour = colors[3],shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 4 + startx + 3*mindist, y = lml_ISR_enzyme$lml[5], colour = colors[3],shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 0 + startx + 4*mindist, y = lml_BSR_enzyme$lml[1], colour = colors[4],shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 1 + startx + 4*mindist, y = lml_BSR_enzyme$lml[2], colour = colors[4],shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 2 + startx + 4*mindist, y = lml_BSR_enzyme$lml[3], colour = colors[4],shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 3 + startx + 4*mindist, y = lml_BSR_enzyme$lml[4], colour = colors[4],shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 4 + startx + 4*mindist, y = lml_BSR_enzyme$lml[5], colour = colors[4],shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 0 + startx + 5*mindist, y = lml_RIF_enzyme$lml[1], colour = colors[5],shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 1 + startx + 5*mindist, y = lml_RIF_enzyme$lml[2], colour = colors[5],shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 2 + startx + 5*mindist, y = lml_RIF_enzyme$lml[3], colour = colors[5],shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 3 + startx + 5*mindist, y = lml_RIF_enzyme$lml[4], colour = colors[5],shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 4 + startx + 5*mindist, y = lml_RIF_enzyme$lml[5], colour = colors[5],shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 0 + startx + 6*mindist, y = lml_RIR_enzyme$lml[1], colour = colors[6],shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 1 + startx + 6*mindist, y = lml_RIR_enzyme$lml[2], colour = colors[6],shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 2 + startx + 6*mindist, y = lml_RIR_enzyme$lml[3], colour = colors[6],shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 3 + startx + 6*mindist, y = lml_RIR_enzyme$lml[4], colour = colors[6],shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 4 + startx + 6*mindist, y = lml_RIR_enzyme$lml[5], colour = colors[6],shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 0 + startx + 2*mindist, y = lml_thames_enzyme$lml[1], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  annotate("point", x = 1 + startx + 2*mindist, y = lml_thames_enzyme$lml[2], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  annotate("point", x = 2 + startx + 2*mindist, y = lml_thames_enzyme$lml[3], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  annotate("point", x = 3 + startx + 2*mindist, y = lml_thames_enzyme$lml[4], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  annotate("point", x = 4 + startx + 2*mindist, y = lml_thames_enzyme$lml[5], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  geom_vline(xintercept=2.5) +
  geom_vline(xintercept=3.5) +
  geom_vline(xintercept=4.5) +
  geom_vline(xintercept=5.5) +
  geom_vline(xintercept=6.5) + labs(title="enzyme dataset",y=" ")+coord_cartesian(ylim=c(-90,-70),xlim = c(1.7,6.3))

annotations_enzyme <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("B"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1)) 

plotenzyme = plotenzyme + theme(axis.title=element_text(size=10)) +
  geom_text(data=annotations_enzyme,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),size = 5)


plotacidity = ggplot(df_acidity[(df_acidity$sampler=="JAGS")&(df_acidity$relabalg=="ECR")&(df_acidity$ellipsealg=="standard")&(df_acidity$thamesalg=="permutations")&(df_acidity$estimate=="estim"),],
                    aes(x=G,y=value)) + annotate("point", x = 0 + startx + 0*mindist, y = lml_true_acidity$lml[1], colour = colors[1],shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 1 + startx + 0*mindist, y = lml_true_acidity$lml[2], colour = colors[1],shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 2 + startx + 0*mindist, y = lml_true_acidity$lml[3], colour = colors[1],shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 3 + startx + 0*mindist, y = lml_true_acidity$lml[4], colour = colors[1],shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 4 + startx + 0*mindist, y = lml_true_acidity$lml[5], colour = colors[1],shape=shapes[1], size=size_ggplot) +
  annotate("point", x = 0 + startx + 1*mindist, y = lml_ISF_acidity$lml[1], colour = colors[2],shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 1 + startx + 1*mindist, y = lml_ISF_acidity$lml[2], colour = colors[2],shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 2 + startx + 1*mindist, y = lml_ISF_acidity$lml[3], colour = colors[2],shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 3 + startx + 1*mindist, y = lml_ISF_acidity$lml[4], colour = colors[2],shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 4 + startx + 1*mindist, y = lml_ISF_acidity$lml[5], colour = colors[2],shape=shapes[2], size=size_ggplot) +
  annotate("point", x = 0 + startx + 3*mindist, y = lml_ISR_acidity$lml[1], colour = colors[3],shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 1 + startx + 3*mindist, y = lml_ISR_acidity$lml[2], colour = colors[3],shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 2 + startx + 3*mindist, y = lml_ISR_acidity$lml[3], colour = colors[3],shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 3 + startx + 3*mindist, y = lml_ISR_acidity$lml[4], colour = colors[3],shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 4 + startx + 3*mindist, y = lml_ISR_acidity$lml[5], colour = colors[3],shape=shapes[3], size=size_ggplot) +
  annotate("point", x = 0 + startx + 4*mindist, y = lml_BSR_acidity$lml[1], colour = colors[4],shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 1 + startx + 4*mindist, y = lml_BSR_acidity$lml[2], colour = colors[4],shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 2 + startx + 4*mindist, y = lml_BSR_acidity$lml[3], colour = colors[4],shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 3 + startx + 4*mindist, y = lml_BSR_acidity$lml[4], colour = colors[4],shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 4 + startx + 4*mindist, y = lml_BSR_acidity$lml[5], colour = colors[4],shape=shapes[4], size=size_ggplot) +
  annotate("point", x = 0 + startx + 5*mindist, y = lml_RIF_acidity$lml[1], colour = colors[5],shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 1 + startx + 5*mindist, y = lml_RIF_acidity$lml[2], colour = colors[5],shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 2 + startx + 5*mindist, y = lml_RIF_acidity$lml[3], colour = colors[5],shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 3 + startx + 5*mindist, y = lml_RIF_acidity$lml[4], colour = colors[5],shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 4 + startx + 5*mindist, y = lml_RIF_acidity$lml[5], colour = colors[5],shape=shapes[5], size=size_ggplot) +
  annotate("point", x = 0 + startx + 6*mindist, y = lml_RIR_acidity$lml[1], colour = colors[6],shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 1 + startx + 6*mindist, y = lml_RIR_acidity$lml[2], colour = colors[6],shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 2 + startx + 6*mindist, y = lml_RIR_acidity$lml[3], colour = colors[6],shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 3 + startx + 6*mindist, y = lml_RIR_acidity$lml[4], colour = colors[6],shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 4 + startx + 6*mindist, y = lml_RIR_acidity$lml[5], colour = colors[6],shape=shapes[6], size=size_ggplot) +
  annotate("point", x = 0 + startx + 2*mindist, y = lml_thames_acidity$lml[1], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  annotate("point", x = 1 + startx + 2*mindist, y = lml_thames_acidity$lml[2], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  annotate("point", x = 2 + startx + 2*mindist, y = lml_thames_acidity$lml[3], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  annotate("point", x = 3 + startx + 2*mindist, y = lml_thames_acidity$lml[4], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  annotate("point", x = 4 + startx + 2*mindist, y = lml_thames_acidity$lml[5], colour = colors[7], shape=shapes[7], size=size_ggplot) +
  geom_vline(xintercept=2.5) +
  geom_vline(xintercept=3.5) +
  geom_vline(xintercept=4.5) +
  geom_vline(xintercept=5.5) +
  geom_vline(xintercept=6.5) + labs(title="acidity dataset",y=" ") + coord_cartesian(ylim=c(-215,-195),xlim = c(1.7,6.3))


annotations_acidity <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("C"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1)) 
plotacidity = plotacidity + theme(axis.title=element_text(size=10)) +
  geom_text(data=annotations_acidity,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),size = 5)


grid.arrange(plotgalaxies,plotenzyme,plotacidity,nrow=3)

# Create a separate plot to act as the legend
data <- data.frame(
  x = c(1, 2, 3, 4, 5, 6, 7),
  y = c(3, 4, 2, 5, 5, 6, 6),
  groupp = c("A", "B", "C", "D","E","F","G")
)


legend_plot <- ggplot(data, aes(x = x, y = y, color = group)) +
  geom_point(data=data[1,],aes(color="BSF"),size = size_ggplot,shape=shapes[1]) + 
  geom_point(data=data[2,],aes(color="ISF"),size = size_ggplot,shape=shapes[2]) +
  geom_point(data=data[2,],aes(color="ISR"),size = size_ggplot,shape=shapes[3]) +
  geom_point(data=data[2,],aes(color="BSR"),size = size_ggplot,shape=shapes[4]) +
  geom_point(data=data[2,],aes(color="RIF"),size = size_ggplot,shape=shapes[5]) +
  geom_point(data=data[2,],aes(color="RIR"),size = size_ggplot,shape=shapes[6]) +
  geom_point(data=data[2,],aes(color="THAMES"),size = size_ggplot,shape=shapes[7]) + 
  scale_color_manual(values = c("THAMES" = colors[7],"BSF" = colors[1], 
                                "ISF" = colors[2],
                                "ISR" = colors[3],
                                "BSR" = colors[4],
                                "RIF" = colors[5],
                                "RIR" = colors[6])) +
  theme(
    legend.position = "bottom", # Legend positioning
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
  ) + labs(color = "estimators")


library(cowplot)

legend=get_plot_component(legend_plot,"guide-box-bottom")
final_plot = plot_grid(plotgalaxies,plotenzyme,plotacidity,nrow=1)
ggsave(filename="atelier/thames_richardson_green.pdf",final_plot,height=2.5,width=8)
ggsave(filename="atelier/thames_richardson_green_legend.pdf",legend,height=2.5,width=8)

df_galaxies_1sim = read.csv('data/res_galaxies_1sim.csv')
load(paste0('data/res_galaxies_1sim_graphs','.Rda'))

#plot(2:15,df_galaxies_1sim[df_galaxies_1sim$estimate=="estim",]$value,xlab="G",ylab="log marginal likelihood")
1+which.max(df_galaxies_1sim[df_galaxies_1sim$estimate=="estim",]$value)

library(igraph)
co = sapply(1:14, function(s) sum(V(graphs[[s]])$color=="blue")-((s+1)-sum(V(graphs[[s]])$color=="blue")))
1+which.max(co)
#plot(2:15,co,xlab="G",ylab="CO")
rbind(2:15,df_galaxies_1sim[df_galaxies_1sim$estimate=="estim",]$value,co)

# proportion of permlen
df_galaxies_1sim[df_galaxies_1sim$estimate=="permlen",]$value/factorial(2:15)

igraphs = list()
#library(ig2gg)
par(mfrow=c(2,7))

for(i in 1:14){
  par(mai=c(0,0,0,0))
  V(graphs[[i]])$color = rep("white",length(V(graphs[[i]])$color))
  V(graphs[[i]])$label = rep(" ",length(V(graphs[[i]])$color))
  
  plot(graphs[[i]],xlim = c(-1,1))
  par(mai=c(0.4,0.4,0.2,0))
  title(paste0("G = ",i+1),adj = 0, line = -0,cex.main = 2)
  text(-1,1,LETTERS[i],xpd=TRUE,cex=1.5)
  abline(v=1.1)
}

graph_example_01 = graph_from_adjacency_matrix(diag(2))

# library(igraph)
# g <- make_ring(1, directed = TRUE)
# g[1]=0
# al <- as_adj_list(g, mode = "out")
# g2 <- graph_from_adj_list(al)
# V(g2)$color="white"
# V(g2)$label=" "
# 
# pdf("atelier/galaxies_overlapgraph_G1.pdf")
# plot(g2)
# dev.off()
# 
# which.max(results_galaxies_1sim$df_output[results_galaxies_1sim$df_output$estimate=="estim",]$value)

