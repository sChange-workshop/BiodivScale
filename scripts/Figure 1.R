### HEADER #####################################################################
##' @title Figure 1 of Chase et al. Species richness change across spatial scales. Oikos
##'
##' @author Patrick Thompson
##' @contact patrick.thompson@zoology.ubc.ca
##' 
##' @date 2018
##' 
##' @description Creates figure 1 from manuscript
##' @note run non_equilibrium_model.r and neutral_model.R to generate the data first
##' 
##' 
##' ################################################################################

library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggExtra)

load("./workspace/non_equil.RData")
load("./workspace/neutral_model_out.RData")

B<-LRR %>% 
  mutate(scale = replace(scale, scale ==2, 10),
         scale = replace(scale, scale == 3, 100),
         scale = replace(scale, scale == 4, 300)) %>% 
  filter(rep == 1) %>% 
  mutate(dispersal_change = replace(dispersal_change, dispersal_change==0, "= dispersal")) %>% 
  mutate(dispersal_change = replace(dispersal_change, dispersal_change=="-0.2", "+ dispersal")) %>% 
  mutate(dispersal_change = replace(dispersal_change, dispersal_change=="0.05", "- dispersal")) %>% 
  mutate(dispersal_change = factor(dispersal_change, levels = c("- dispersal","= dispersal", "+ dispersal"), ordered = TRUE)) %>% 
  mutate(abundance_change = replace(abundance_change, abundance_change==0, "= carrying capacity")) %>% 
  mutate(abundance_change = replace(abundance_change, abundance_change=="-10", "+ carrying capacity")) %>%
  mutate(abundance_change = replace(abundance_change, abundance_change=="10", "- carrying capacity")) %>%
  mutate(abundance_change = factor(abundance_change, levels = rev(c("- carrying capacity", "= carrying capacity", "+ carrying capacity")), ordered = TRUE)) %>% 
  ggplot(aes(x=factor(scale), y = LRR, color = abundance_change, group = abundance_change))+
  geom_hline(yintercept = 0,linetype=2)+
  facet_wrap(~dispersal_change)+
  geom_boxplot(aes(group = interaction(abundance_change, scale)))+
  stat_summary(fun.y=median, position=position_dodge(width=0.75), geom="line", aes(group=abundance_change))  + 
  scale_color_manual(values = c("red", "grey30","dodgerblue"), name = "")+
  theme_bw()+
  removeGrid()+
  #scale_x_log10(breaks = c(1,10,100,300))+
  xlab("Spatial scale (# of communities)")+
  ylab("Species richness change (LRS)")+
  theme(legend.justification=c(1.05,-0.05), legend.position=c(1,0))+
  theme(legend.title=element_blank())+
  theme(strip.background = element_rect(colour=NA, fill=NA))


A<-ggplot(filter(Change.df, time == 100),aes(x=factor(Scale_size),y=log_ratio,color=as.factor(Local_loss),group=interaction(Scale_size,Local_loss)))+
  geom_hline(yintercept = 0,linetype=2)+
  #geom_line(data = medians, position=position_dodge(width=0.25),aes(x= as.factor(Scale_size), y=log_ratio,color=as.factor(Scenario),group=Scenario))+
  stat_summary(fun.y=median, position=position_dodge(width=0.75), geom="line", aes(group=Local_loss))  + 
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~Global_spread_text)+
  scale_color_manual(values = rev(brewer.pal(n = 11,name = "Spectral")), name = "local loss")+
  removeGrid()+
  xlab("Spatial scale (# of communities)")+
  #scale_x_log10(breaks = c(1,10,100,1000))+
  ylab("Species richness change (LRS)")+
  theme(legend.justification=c(1.05,1.05), legend.position=c(1,1))+
  theme(legend.title=element_blank())+
  theme(strip.background = element_rect(colour=NA, fill=NA))

library(cowplot)

plot_grid(A,B,nrow = 2, labels = c("(A)", "(B)"))
ggsave("./figures/Figure 1.pdf",width = 12, height =8)