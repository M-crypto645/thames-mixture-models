# simulate from the Gaussian mixture model with known proportions and variances
# and compute the marginal likelihood estimators
rm(list=ls())
source("functions/galaxies_funcs_squares.R")
source("functions/true_marglik_funcs.R") 
source("functions/pipeline.R")
library(reshape2)
library(combinat)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(gridExtra)

Sigmahatinv = cbind(c(0.9,-0.37),c(-0.37,1.74))
thetahat = c(1.65,4.35)
c <- 2+1                              # radius

# transform for visualization
eig <- eigen(Sigmahatinv)
Q <- eig$vectors
D <- diag(eig$values)
angle <- seq(0, 2*pi, length.out = 100)
circle <- rbind(cos(angle), sin(angle))
Dinv_sqrt <- diag(1 / sqrt(eig$values))
ellipse <- Q %*% Dinv_sqrt %*% circle * sqrt(c)
ellipse[1, ] <- ellipse[1, ] + thetahat[1]
ellipse[2, ] <- ellipse[2, ] + thetahat[2]

# Step 5: Make data frame for ggplot
df_ellipse <- data.frame(x = ellipse[1, ], y = ellipse[2, ])

#mutate(V3=case_when(V1>V2~"1",T~"0"))
df_ellipse$V3=0

df = as.data.frame(rbind(c(1,4),
                         c(3,4.3),
                         c(2,4.5),
                         c(1.3,4.2),
                         c(0.9,4.1),
                         c(3.9,1.1),
                         c(4.5,2.5),
                         c(4,1.5),
                         c(5.1,1),
                         c(6,1))) %>% 
  mutate(V3=paste0(1:10))

plot.title.size = 23
axis.title.x.size = 32
axis.title.y.size = 32
axis.text.size = 28
label.size = 10
plot_psi2 =  df %>% 
  mutate(V4=case_when(V1>V2~V1,T~V2),V5=case_when(V2>V1~V1,T~V2)) %>% 
  ggplot(aes(x=V4,y=V5,label=V3)) +
  geom_text_repel(aes(segment.colour = "grey50"), size = label.size) +
  geom_point() +
  geom_abline(slope=1,intercept=0,color="red") +
  geom_path(data=df_ellipse,aes(x,y),color = "brown", size = 1.2)+
  xlim(-3,8) +
  scale_color_discrete(label=c(expression(mu[1]),expression(mu[1])))+
  ylim(-3,8) +
  labs(x=expression(mu[1]),y=expression(mu[2]),
       title=expression(psi[2]~"(0 points in the ellipsoid)")) +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )

plot_psi1 =  df %>% 
  mutate(V4=case_when(V1>V2~V2,T~V1),V5=case_when(V2>V1~V2,T~V1)) %>% 
  ggplot(aes(x=V4,y=V5,label=V3)) +
  geom_path(data=df_ellipse,aes(x,y),color = "brown", size = 1.2)+
  geom_text_repel(aes(segment.colour = "grey50"), size = label.size) + 
  geom_point() +
  geom_abline(slope=1,intercept=0,color="red") +
  xlim(-3,8) +
  scale_color_discrete(label=c(expression(mu[1]),expression(mu[1])))+
  ylim(-3,8) +
  labs(x=expression(mu[1]),y=expression(mu[2]),
       title=expression(psi[1]~"(9 points in the ellipsoid)")) +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )

plot_P1 =  df %>% 
  mutate(V4=V1,V5=V2) %>% 
  ggplot(aes(x=V4,y=V5,label=V3)) +
  geom_path(data=df_ellipse,aes(x,y),color = "brown", size = 1.2)+
  geom_text_repel(aes(segment.colour = "grey50"), size = label.size) + 
  geom_point() +
  geom_abline(slope=1,intercept=0,color="red") +
  xlim(-3,8) +
  scale_color_discrete(label=c(expression(mu[1]),expression(mu[1])))+
  ylim(-3,8) +
  labs(x=expression(mu[1]),y=expression(mu[2]),
       title=expression(P[1]~"(5 points in the ellipsoid)")) +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )

plot_P2 =  df %>% 
  mutate(V4=V2,V5=V1) %>% 
  ggplot(aes(x=V4,y=V5,label=V3)) +
  geom_path(data=df_ellipse,aes(x,y),color = "brown", size = 1.2)+
  geom_text_repel(aes(segment.colour = "grey50"), size = label.size) +
  geom_point() +
  geom_abline(slope=1,intercept=0,color="red") +
  xlim(-3,8) +
  scale_color_discrete(label=c(expression(mu[1]),expression(mu[1])))+
  ylim(-3,8) +
  labs(x=expression(mu[1]),y=expression(mu[2]),
       title=expression(P[2]~"(4 points in the ellipsoid)")) +
  theme(
    plot.title = element_text(size = plot.title.size),
    axis.title.x = element_text(size = axis.title.x.size),
    axis.title.y = element_text(size = axis.title.y.size),
    axis.text = element_text(size = axis.text.size)
  )



annotations_P1 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("A"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))
annotations_P2 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("B"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))
annotations_Psi1 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("C"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))
annotations_Psi2 <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("D"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1))

letter_size = 10
plot_final = grid.arrange(plot_P1 + geom_text(data=annotations_P1,
                                              aes(x=xpos,
                                                  y=ypos,
                                                  hjust=hjustvar,
                                                  vjust=vjustvar,
                                                  label=annotateText),
                                              size = letter_size),
                          plot_P2 + geom_text(data=annotations_P2,
                                              aes(x=xpos,
                                                  y=ypos,
                                                  hjust=hjustvar,
                                                  vjust=vjustvar,
                                                  label=annotateText),
                                              size = letter_size),
                          plot_psi1  + geom_text(data=annotations_Psi1,
                                                 aes(x=xpos,
                                                     y=ypos,
                                                     hjust=hjustvar,
                                                     vjust=vjustvar,
                                                     label=annotateText),
                                                 size = letter_size),
                          plot_psi2 + geom_text(data=annotations_Psi2,
                                                aes(x=xpos,
                                                    y=ypos,
                                                    hjust=hjustvar,
                                                    vjust=vjustvar,
                                                    label=annotateText),
                                                size = letter_size),
                          ncol=4)
ggsave("atelier/dummy_example.pdf",plot_final,height=5,width=20)
