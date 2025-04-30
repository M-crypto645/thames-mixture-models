# Fit the THAMES to the galaxies, enzyme, and acidity datasets
rm(list=ls())
# if(strsplit(getwd(),"/")[[1]][length(strsplit(getwd(),"/")[[1]])]!="thames_mixtures"){
#   setwd("thames_mixtures")
# }
#source('galaxies_funcs_squares.R')
source('functions/galaxies_funcs_squares.R')
source('functions/pipeline.R')
source("functions/theme.R")
#source('functions/thames_gmm.R')
#source("functions/thames_gmm_funcs.R")
pacman::p_load(rstan,label.switching,combinat)
#options(mc.cores = parallel::detectCores())
# remove.packages("dplyr")
# install.packages("dplyr",dependencies=TRUE)
library(multimode)
library(ggplot2)
library(gridExtra)
library(grid)

G_list = 2:6
dfs_galaxies = list()
for(g in seq_along(G_list)){
  dfs_galaxies[[g]] = read.csv(file=paste0("data/thetaW_galaxies_G",G_list[g],".csv"))
}
names(dfs_galaxies) = paste0(rep("G",6-2+1),2:6)

# floor to indicate cluster assignment
dfs_galaxies$G4$WGrid1 = floor(dfs_galaxies$G4$WGrid1)
dfs_galaxies$G4$WGrid2 = floor(dfs_galaxies$G4$WGrid2)
dfs_galaxies$G4$WGrid3 = floor(dfs_galaxies$G4$WGrid3)
dfs_galaxies$G4$WGrid4 = floor(dfs_galaxies$G4$WGrid4)

plot.title.size = 20
axis.title.x.size = 16
axis.title.y.size = 16
axis.text.size = 14

plot_comp1_G4 = ggplot(arrange(dfs_galaxies$G4,W1),aes(x=mu1,y=sigmasqu1,color=W1)) +
  geom_point(aes(alpha=1)) +
  labs(x=expression(mu[1]),y=expression(sigma[1]^2)) +
  scale_color_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) +
  geom_raster(aes(x = muGrid1, y = sigmasquGrid1, fill = WGrid1),alpha=0.3)+
  scale_fill_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  )+ theme(legend.position = "none") +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )
plot_comp2_G4 = ggplot(arrange(dfs_galaxies$G4,W2),aes(x=mu2,y=sigmasqu2,color=W2)) +
  geom_point(alpha=1) +
  labs(x=expression(mu[2]),y=expression(sigma[2]^2)) +
  scale_color_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) +
  geom_raster(aes(x = muGrid2, y = sigmasquGrid2, fill = WGrid2),alpha=0.18)+
  scale_fill_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) + theme(legend.position = "none") +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )
plot_comp3_G4 = ggplot(arrange(dfs_galaxies$G4,W3),aes(x=mu3,y=sigmasqu3,color=W3)) +
  geom_point(alpha=1) +
  labs(x=expression(mu[3]),y=expression(sigma[3]^2)) +
  scale_color_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) +
  geom_raster(aes(x = muGrid3, y = sigmasquGrid3, fill = WGrid3),alpha=0.18)+
  scale_fill_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) + theme(legend.position = "none") +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )

plot_comp4_G4 = ggplot(arrange(dfs_galaxies$G4,W4),aes(x=mu4,y=sigmasqu4,color=W4)) +
  geom_point(alpha=1) +
  labs(x=expression(mu[4]),y=expression(sigma[4]^2)) +
  scale_color_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) +
  #geom_point(aes(x = muGrid4, y = sigmasquGrid4,color=WGrid4))
  geom_raster(aes(x = muGrid4, y = sigmasquGrid4, fill = WGrid4),alpha=0.18,
            , width = 1.25, height = 1.25)+
  scale_fill_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) + theme(legend.position = "none") +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )

annotations_comp1G4 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("F"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))
annotations_comp2G4 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("G"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))
annotations_comp3G4 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("H"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))
annotations_comp4G4 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c(" I"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))

plot_G4 = grid.arrange(plot_comp1_G4 + geom_text(data=annotations_comp1G4,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       plot_comp2_G4 + geom_text(data=annotations_comp2G4,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       plot_comp3_G4 + geom_text(data=annotations_comp3G4,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       plot_comp4_G4 + geom_text(data=annotations_comp4G4,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       ncol=2,
                       top = textGrob("G = 4", gp = gpar(fontsize = 20, fontface = "bold")))


plot_G4a = grid.arrange(plot_comp1_G4 + geom_text(data=annotations_comp1G4,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       plot_comp2_G4 + geom_text(data=annotations_comp2G4,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       ncol=2,
                       top = textGrob("G = 4", gp = gpar(fontsize = 20, fontface = "bold")))

plot_G4b = grid.arrange(plot_comp3_G4 + geom_text(data=annotations_comp3G4,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       plot_comp4_G4 + geom_text(data=annotations_comp4G4,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       ncol=2)

# ggsave("atelier/vis_galaxy_G4.pdf",plot_G4,height=2.5,width=8)

# floor to indicate cluster assignment
dfs_galaxies$G3$WGrid1 = floor(dfs_galaxies$G3$WGrid1)
dfs_galaxies$G3$WGrid2 = floor(dfs_galaxies$G3$WGrid2)
dfs_galaxies$G3$WGrid3 = floor(dfs_galaxies$G3$WGrid3)

# plot.title.size = 20
# axis.title.x.size = 16
# axis.title.y.size = 16
# axis.text.size = 14

plot_comp1_G3 = ggplot(arrange(dfs_galaxies$G3,W1),aes(x=mu1,y=sigmasqu1,color=W1)) +
  geom_point(aes(alpha=1)) +
  labs(x=expression(mu[1]),y=expression(sigma[1]^2)) +
  scale_color_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) +
  geom_raster(aes(x = muGrid1, y = sigmasquGrid1, fill = WGrid1),alpha=0.3)+
  scale_fill_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) + theme(legend.position = "none") +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )
plot_comp2_G3 = ggplot(arrange(dfs_galaxies$G3,W2),aes(x=mu2,y=sigmasqu2,color=W2)) +
  geom_point(alpha=1) +
  labs(x=expression(mu[2]),y=expression(sigma[2]^2)) +
  scale_color_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) +
  geom_raster(aes(x = muGrid2, y = sigmasquGrid2, fill = WGrid2),alpha=0.15)+
  scale_fill_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) + theme(legend.position = "none") +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )
plot_comp3_G3 = ggplot(arrange(dfs_galaxies$G3,W3),aes(x=mu3,y=sigmasqu3,color=W3)) +
  geom_point(alpha=1) +
  labs(x=expression(mu[3]),y=expression(sigma[3]^2)) +
  scale_color_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) +
  geom_raster(aes(x = muGrid3, y = sigmasquGrid3, fill = WGrid3),alpha=0.13)+
  scale_fill_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) + theme(legend.position = "none") +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )

annotations_comp1G3 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("C"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))
annotations_comp2G3 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("D"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))
annotations_comp3G3 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("E"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))

plot_G3 = grid.arrange(plot_comp1_G3 + geom_text(data=annotations_comp1G3,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       plot_comp2_G3 + geom_text(data=annotations_comp2G3,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       plot_comp3_G3 + geom_text(data=annotations_comp3G3,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       ncol=3, 
                       top = textGrob("G = 3", gp = gpar(fontsize = 20, fontface = "bold")))
# ggsave("atelier/vis_galaxy_G3.pdf",plot_G3,height=2.5,width=8)

# floor to indicate cluster assignment
dfs_galaxies$G2$WGrid1 = floor(dfs_galaxies$G2$WGrid1)
dfs_galaxies$G2$WGrid2 = floor(dfs_galaxies$G2$WGrid2)

plot_comp1_G2 = ggplot(arrange(dfs_galaxies$G2,W1),aes(x=mu1,y=sigmasqu1,color=W1)) +
  geom_point(aes(alpha=1)) +
  labs(x=expression(mu[1]),y=expression(sigma[1]^2)) +
  scale_color_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) +
  geom_raster(aes(x = muGrid1, y = sigmasquGrid1, fill = WGrid1),alpha=0.3)+
  scale_fill_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) + theme(legend.position = "none") +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )
plot_comp2_G2 = ggplot(arrange(dfs_galaxies$G2,W2),aes(x=mu2,y=sigmasqu2,color=W2)) +
  geom_point(alpha=1) +
  labs(x=expression(mu[2]),y=expression(sigma[2]^2)) +
  scale_color_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) +
  geom_raster(aes(x = muGrid2, y = sigmasquGrid2, fill = WGrid2),alpha=0.15)+
  scale_fill_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) + theme(legend.position = "none") +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )

annotations_comp1G2 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("A"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))
annotations_comp2G2 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("B"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))

library(grid)
plot_G2 = grid.arrange(plot_comp1_G2 + geom_text(data=annotations_comp1G2,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       plot_comp2_G2 + geom_text(data=annotations_comp2G2,aes(x=xpos,
                                                                              y=ypos,
                                                                              hjust=hjustvar,
                                                                              vjust=vjustvar,
                                                                              label=annotateText,
                                                                              fill=NULL,
                                                                              color=NULL),size=5),
                       ncol=2, 
                       top = textGrob("G = 2", gp = gpar(fontsize = 20, fontface = "bold")))
# ggsave("atelier/vis_galaxy_G2.pdf",plot_G2,height=2.5,width=8)

#ggsave("atelier/vis_galaxy_full.pdf",
#       grid.arrange(plot_G2,plot_G3,plot_G4,nrow=3))
ggsave("atelier/vis_galaxy_full.pdf",
       grid.arrange(plot_G2,plot_G3,
                    plot_G4a,plot_G4b,nrow=4))

library(cowplot)
plot_legend = ggplot(arrange(dfs_galaxies$G4,W1),aes(x=mu1,y=sigmasqu1,color=W1)) +
  geom_point() +
  labs(x=expression(mu[1]),y=expression(sigma[1]^2)) +
  scale_color_gradient2(
    low = "cyan",
    mid = "black",
    high = "red",
    midpoint = 2,
    limits = c(1,3.5)
  ) + theme(legend.position = "bottom",
            legend.key.width = unit(3, "cm")) + labs(color="W\n ")
ggsave(filename="atelier/vis_galaxy_legend.pdf",
       cowplot::get_plot_component(plot_legend,"guide-box-bottom"),height=2.5,width=8)
