### HEADER #####################################################################
##' @title Figure 2 of Chase et al. Species richness change across spatial scales. Oikos https://doi.org/10.1111/oik.05968
##' 
##' @date 2018
##' 
##' @description Creates figure 2 from manuscript
##' 
##' 
##' ################################################################################

library(tidyverse)
library(ggExtra)
library(cowplot)
library(scales)

#Panama corals####
corals <- read.csv("./data_for_dryad/panama_corals.csv")

corals$Type <- "Central American corals"

coral_inset<-corals %>%
  gather(key = year, value = S, S1:S2) %>% 
  ggplot(aes(x = scale, y = S, group = year, linetype = factor(year)))+
  geom_smooth(method = "nls", formula = y ~ a * x^b, color = 1, se = F, size = 0.5)+
  scale_linetype(name = "", guide = F)+
  xlab(expression(paste("Scale (",m^2,")",sep="")))+
  ylab("Richness change per year")+
  ylab("Richness")+
  theme(legend.justification=c(1,0), legend.position=c(1,0), text = element_text(size = 8))+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9))


coral_fig<-ggplot(corals,aes(x=scale,y=S_change_per_year))+
  geom_hline(yintercept = 0, linetype=1)+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = 1, se = FALSE)+
  stat_quantile(formula = y ~ poly(x, 2),quantiles = c(0.1,0.9), color="black", linetype=2)+
  theme_bw()+
  facet_wrap(~Type)+
  theme(strip.background = element_rect(colour="NA", fill="NA"),strip.text.x = element_text(size = 12))+
  xlab(expression(paste("Scale (",m^2,")",sep="")))+
  ylab("Richness change per year")+
  scale_x_log10(breaks = c(seq(1,9,by=1),seq(10,90,by = 10),seq(100,900,by = 100)),
                labels = c(1,rep("",8),10,rep("", 8),100,rep("", 8)))+
  removeGrid()

Fig2a<- ggdraw()+
  draw_plot(coral_fig)+
  draw_plot(coral_inset, x = 0.66, y = 0.58, width = 0.33, height = 0.35)

#EU plants####
European_plants_LRR <- read.csv("./data_for_dryad/EU_plants.csv")
European_plants_LRR$Type <- "European plants"

plant_inset <- European_plants_LRR %>% 
  ungroup() %>% 
  select(Area, S1, S2) %>% 
  gather(key = time, value = S, S1:S2) %>% 
  ggplot(aes(x = Area, y = S, group = time, linetype = factor(time)))+
  geom_smooth(method = "nls", formula = y ~ a * x^b, color = 1, se = F, size = 0.5)+
  scale_linetype(name = "", guide = F)+
  xlab(expression(paste("Scale (",km^2,")",sep="")))+
  ylab("Richness change per year")+
  ylab("Richness")+
  theme(legend.justification=c(1,0), legend.position=c(1,0), text = element_text(size = 8))+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9))+
  scale_x_continuous(breaks = seq(0,5e6, by = 2e6))

plant_fig<-ggplot(European_plants_LRR,aes(x = Area,y=LRR))+
  geom_hline(yintercept = 0, linetype=1)+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = 1, se = FALSE)+
  stat_quantile(formula = y ~ poly(x, 2),quantiles = c(0.1,0.9), color="black", linetype=2)+
  theme_bw()+
  facet_wrap(~Type)+
  theme(strip.background = element_rect(colour="NA", fill="NA"),strip.text.x = element_text(size = 12))+
  scale_x_log10(breaks = c(seq(1e4,9e4,by=1e4),seq(1e5,9e5,by = 1e5),seq(1e6,9e6,by = 1e6)),
                labels = c(1e4,rep("",8),1e5,rep("", 8),1e6,rep("", 8)))+
  xlab(expression(paste("Scale (",km^2,")",sep="")))+
  ylab("LRS")+
  removeGrid()

Fig2d<- ggdraw()+
  draw_plot(plant_fig)+
  draw_plot(plant_inset, x = 0.66, y = 0.58, width = 0.33, height = 0.35)

#Hawaiian birds####

hawaii <- read.csv("./data_for_dryad/hawaii_birds.csv")
hawaii$Type <- "Hawaiian birds"
hawaii$LRR<-log(hawaii$SpeciesNow/hawaii$SpeciesBefore)

hawaii_inset <- hawaii %>% 
  ungroup() %>% 
  select(Area_km_2, SpeciesBefore, SpeciesNow) %>% 
  gather(key = time, value = S, SpeciesBefore:SpeciesNow) %>% 
  ggplot(aes(x = Area_km_2, y = S, group = time, linetype = factor(time)))+
  geom_smooth(method = "nls", formula = y ~ a * x^b, color = 1, se = F, size = 0.5)+
  scale_linetype(name = "", guide = F)+
  xlab(expression(paste("Scale (",km^2,")",sep="")))+
  ylab("Richness change per year")+
  ylab("Richness")+
  theme(legend.justification=c(1,0), legend.position=c(1,0), text = element_text(size = 8))+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9))

hawaii_fig<-ggplot(hawaii,aes(x = Area_km_2, y=LRR))+
  geom_hline(yintercept = 0, linetype=1)+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = 1, se = FALSE)+
  stat_quantile(formula = y ~ poly(x, 2),quantiles = c(0.1,0.9), color="black", linetype=2)+
  theme_bw()+
  facet_wrap(~Type)+
  theme(strip.background = element_rect(colour="NA", fill="NA"),strip.text.x = element_text(size = 12))+
  ylab("LRS")+
  xlab(expression(paste("Scale (",km^2,")",sep="")))+
  scale_x_log10(breaks = c(seq(100,900,by=100),seq(1000,10000,by = 1000)),
                labels = c(100,rep("",8),1000,rep("", 8),10000))+
  removeGrid()

Fig2c<- ggdraw()+
  draw_plot(hawaii_fig)+
  draw_plot(hawaii_inset, x = 0.66, y = 0.58, width = 0.33, height = 0.35)

#BBS####
bbs <- read.csv("./data_for_dryad/bbs_data.csv")
bbs$Type <- "North American birds"

bbs_inset <- bbs %>% 
  ungroup() %>% 
  select(Aream2, S1, S2) %>% 
  gather(key = time, value = S, S1:S2) %>% 
  ggplot(aes(x = Aream2*10e-6, y = S, group = time, linetype = factor(time)))+
  geom_smooth(method = "nls", formula = y ~ a * x^b, color = 1, se = F, size = 0.5)+
  scale_linetype(name = "", guide = F)+
  scale_x_continuous(breaks = c(0,1e4, 2e4))+
  xlab(expression(paste("Scale (",km^2,")",sep="")))+
  ylab("Richness change per year")+
  ylab("Richness")+
  theme(legend.justification=c(1,0), legend.position=c(1,0), text = element_text(size = 8))+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9))

bbs_fig<-ggplot(bbs,aes(x = Aream2*10e-6,y=slope))+
  geom_hline(yintercept = 0, linetype=1)+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = 1, se = FALSE)+
  stat_quantile(formula = y ~ poly(x, 2),quantiles = c(0.1,0.9), color="black", linetype=2)+
  theme_bw()+
  facet_wrap(~Type)+
  theme(strip.background = element_rect(colour="NA", fill="NA"),strip.text.x = element_text(size = 12))+
  ylab("Richness change per year")+
  xlab(expression(paste("Scale (",km^2,")",sep="")))+
  scale_x_log10(breaks = c(seq(1e2,1e3,by=1e2),seq(1e3,1e4,by = 1e3),seq(1e4,1e5,by = 1e4)),
                labels = c(1e2,rep("",9),1e3,rep("", 9),1e4,rep("", 9)), limits = c(1e2, 2.5e4))+
  removeGrid()

Fig2b<- ggdraw()+
  draw_plot(bbs_fig)+
  draw_plot(bbs_inset, x = 0.66, y = 0.58, width = 0.33, height = 0.35)

plot_grid(Fig2a,Fig2b,Fig2c,Fig2d,nrow = 2,labels = c("(A)","(B)","(C)","(D)"))
ggsave("./figures/Figure 2.pdf",height = 6*1.5, width = 8*1.5)

