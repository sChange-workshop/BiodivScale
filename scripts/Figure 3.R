### HEADER #####################################################################
##' @title Figure 3 of Chase et al. Species richness change across spatial scales. Oikos https://doi.org/10.1111/oik.05968
##'
##' @date 2018
##' 
##' @description Creates figure 3 from manuscript
##' 
##' ################################################################################

library(tidyverse)
library(ggExtra)
library(cowplot)

all.data <- read.csv("./data_for_dryad/all.data.csv")

bird.data<-all.data %>% 
  filter(Scale >0) %>% 
  #filter(!is.na(Scale)) %>% 
  dplyr::filter(BioType == "terrestrial" | BioType == "global") %>% 
  filter(Taxa == "birds")

Fig3a <- ggplot(bird.data, aes(x=Scale, y=log_ratio))+
  facet_wrap(~Taxa)+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_point(pch=21, aes(fill = StudyType))+ 
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = 1, se = FALSE)+
  stat_quantile(formula = y ~ poly(x,2),quantiles = c(0.1,0.9), color="black", linetype=2)+
  scale_fill_grey(name = "")+
  scale_x_log10(limits=c(1e-8,1.5e+8),breaks=c(1e-8,1e-4,1e0,1e4,1e8))+
  scale_y_continuous(limits = c(-1.8,2.1))+
  theme_bw()+
  removeGrid()+
  theme(strip.background = element_rect(colour="NA", fill="NA"),strip.text.x = element_text(size = 12))+
  ylab("LRS")+
  xlab(expression(paste("Scale (",km^2,")",sep="")))+
  theme(legend.justification=c(0,1), legend.position=c(0.01,0.99))


plant.data<-all.data %>% 
  filter(Scale >0) %>% 
  #filter(!is.na(Scale)) %>% 
  dplyr::filter(BioType == "terrestrial" | BioType == "global") %>% 
  filter(Taxa == "plants")

plant.data$Taxa<-"terrestrial plants"

Fig3b <- ggplot(plant.data, aes(x=Scale, y=log_ratio))+
  facet_wrap(~Taxa)+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_point(pch=21, aes(fill = StudyType))+ 
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = 1, se = FALSE)+
  stat_quantile(formula = y ~ poly(x,2),quantiles = c(0.1,0.9), color="black", linetype=2)+
  scale_fill_grey(guide=FALSE)+
  scale_x_log10(limits=c(1e-8,1.5e+8),breaks=c(1e-8,1e-4,1e0,1e4,1e8))+
  scale_y_continuous(limits = c(-1.8,2.1))+
  theme_bw()+
  removeGrid()+
  theme(strip.background = element_rect(colour="NA", fill="NA"),strip.text.x = element_text(size = 12))+
  ylab("LRS")+
  xlab(expression(paste("Scale (",km^2,")",sep="")))

mammal.data<-all.data %>% 
  filter(Scale >0) %>% 
  #filter(!is.na(Scale)) %>% 
  dplyr::filter(BioType == "terrestrial" | BioType == "global") %>% 
  filter(Taxa == "mammals")

mammal.data$Taxa<-"terrestrial mammals"

Fig3c<- ggplot(mammal.data,aes(x=Scale, y=log_ratio))+
  facet_wrap(~Taxa)+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_point(pch=21, aes(fill = StudyType))+ 
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = 1, se = FALSE)+
  stat_quantile(formula = y ~ poly(x,2),quantiles = c(0.1,0.9), color="black", linetype=2)+
  scale_fill_grey(guide=FALSE)+
  scale_x_log10(limits=c(1e-8,1.5e+8),breaks=c(1e-8,1e-4,1e0,1e4,1e8))+
  scale_y_continuous(limits = c(-1.8,2.1))+
  theme_bw()+
  removeGrid()+
  theme(strip.background = element_rect(colour="NA", fill="NA"),strip.text.x = element_text(size = 12))+
  ylab("LRS")+
  xlab(expression(paste("Scale (",km^2,")",sep="")))

fish.data<-all.data %>% 
  filter(Scale >0) %>% 
  #filter(!is.na(Scale)) %>% 
  dplyr::filter(BioType == "marine" | BioType == "global") %>% 
  filter(Taxa == "fish")

fish.data$Taxa<-"marine fish"

Fig3d<-ggplot(fish.data,aes(x=Scale, y=log_ratio))+
  facet_wrap(~Taxa)+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_point(pch=21, aes(fill = StudyType))+ 
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = 1, se = FALSE)+
  stat_quantile(formula = y ~ poly(x,2),quantiles = c(0.1,0.9), color="black", linetype=2)+
  scale_fill_grey(guide=FALSE)+
  scale_x_log10(limits=c(1e-8,1.5e+8),breaks=c(1e-8,1e-4,1e0,1e4,1e8))+
  scale_y_continuous(limits = c(-1.8,2.1))+
  theme_bw()+
  removeGrid()+
  theme(strip.background = element_rect(colour="NA", fill="NA"),strip.text.x = element_text(size = 12))+
  ylab("LRS")+
  xlab(expression(paste("Scale (",km^2,")",sep="")))

plot_grid(Fig3a,Fig3b,Fig3c,Fig3d,nrow = 2,labels = c("(A)","(B)","(C)","(D)"))
ggsave("./figures/Figure 3.pdf",height = 6, width = 8, dpi = 300)
