### HEADER #####################################################################
##' @title Non equilibirum biodiversity change simulation - 4 scales
##'
##' @author Patrick Thompson
##' @contact patrick.thompson@zoology.ubc.ca
##' 
##' @date 26 Feb 2016
##' 
##' @description Simulates scale dependent biodiversity change from local to global scales. 
##' Consists of 4 nested sets of species pools. 
##' Between time 1 and time 2, a certain proportion of species are lost from each local community 
##' and there is a certain chance of a species arriving from the global pool. 
##' These changes scale up to deterime biodiversity change at larger scales.
##' 
##' @note Species loss can either be correlated across local sites or uncorrelated. This is set using correlate_loss.
##' 
##' ################################################################################

library(ggplot2)
library(vegan)
library(dplyr)
library(ggExtra)
library(RColorBrewer)
library(viridis)

options(scipen = 9999)

Global_pool<-1000
Regions<-10
Countries_per_region<-10
Provinces_per_country<-10

correlate_loss<-T #set to true if the chance that a species is lost should be correlated across sites
correlated_spread<-F

time_points<-100

extirp_probV<-c(0.1,0.03,0.01)#c(0.1,0.25,0.5)
invasion_probV<-0.001 #c(0.0001,0.001,0.01)
global_spread<-matrix(c(1,1,1,1,10,100,1,100,1000),3,3)  # add a 4 to give it a 4th colum, and could add whatever set of 3 K values I want; keep 1 as baseline, and add larger numbers. dispersal probabililty is weighted by abundance. the key thing is the relationship to the 1.1 is the chance to get a local species. x and y coordinates bring in assumptions about edge effects, distance between patches, etc. complicating.

loss_prob<-rbeta(Global_pool,2,5)
spread_prob<-runif(Global_pool)

jaccard<-function(pm){
  with(pm, {
    (b+c)/(a+b+c)
  })}


betapart<-function (A, B) 
  list(b = sum(A==T & B==F), c = sum(A==F & B==T), a = sum(A==T & B==T))

betapart2<-function (A, B) 
  list(b = colSums(A==T & B==F), c = colSums(A==F & B==T), a = colSums(A==T & B==T))

Change.df<-data.frame()
Species.df<-data.frame()

for(g in 1:ncol(global_spread)){
  print(g)
  for(ext in 1:length(extirp_probV)){
    overall_extirp_prob<-extirp_probV[ext]
    for(inv in 1:length(invasion_probV)){
      invasion_prob<-invasion_probV[inv]
      
      Regional_pools<-array(NA,dim=c(Global_pool,Regions,time_points+1))
      Country_pools<-array(NA,dim=c(Global_pool,Regions,Countries_per_region,time_points+1))
      Province_pools<-array(NA,dim=c(Global_pool,Regions,Countries_per_region,Provinces_per_country,time_points+1))
      for(i in 1:Regions){
        prob<-700 #round(rnorm(1,mean=700,sd=50))
        Regional_pools[,i,1]<-sample(c(rep(0,prob),rep(1,Global_pool-prob)),Global_pool,replace=T)
        for(j in 1:Countries_per_region){
          Country_pools[,i,j,1]<-sample(c(rep(0,prob),rep(1,Global_pool-prob)),Global_pool,replace=T)*Regional_pools[,i,1]
          for(k in 1:Provinces_per_country){
            Province_pools[,i,j,k,1]<-sample(c(rep(0,prob),rep(1,Global_pool-prob)),Global_pool,replace=T)*Country_pools[,i,j,1]
          }
          Country_pools[,i,j,1]<-rowSums(Province_pools[,i,j,,1])>0
        }
        Regional_pools[,i,1]<-rowSums(Country_pools[,i,,1])>0
      }
      
      distance_array<-Province_pools[,,,,1]
      distance_array[,,,]<-1/global_spread[3,g]
      
      for(time in 1:time_points){
        print(paste("time = ", time))
        for(i in 1:Regions){
          distance_array[,i,,]<-1/global_spread[2,g]
          for(j in 1:Countries_per_region){
            distance_array[,i,j,]<-1/global_spread[1,g]
            for(k in 1:Provinces_per_country){
              distance_array[,i,j,k]<-0
              dispersal_pool<-rowSums(distance_array*Province_pools[,,,,time])
              dispersal_pool<-dispersal_pool/sum(dispersal_pool)
              distance_array[,i,j,k]<-1/global_spread[1,g]
              Province_pools[,i,j,k,time+1]<-Province_pools[,i,j,k,time]
              disp<-rbinom(Global_pool,1,prob=invasion_prob*dispersal_pool/mean(dispersal_pool))>0
              Province_pools[disp,i,j,k,time+1]<-1#spread step
              # if(correlated_spread==T){
              #   Province_pools[rbinom(Global_pool,1,prob=invasion_prob*spread_prob*2)*rowSums(Regional_pools[,,1])>0,i,j,k,2]<-1
              # } else{
              #   Province_pools[rbinom(Global_pool,1,prob=invasion_prob)*rowSums(Regional_pools[,,1])>0,i,j,k,2]<-1}#spread step
              if(correlate_loss==T){
                Province_pools[sample(Global_pool,overall_extirp_prob*Global_pool,prob = loss_prob,replace=F),i,j,k,time+1]<-0 #loss step
              } else {Province_pools[sample(Global_pool,overall_extirp_prob*Global_pool,replace=F),i,j,k,time+1]<-0} #loss step
            }
            Country_pools[,i,j,time+1]<-rowSums(Province_pools[,i,j,,time+1])>0
            distance_array[,i,j,]<-1/global_spread[2,g]
          }
          Regional_pools[,i,time+1]<-rowSums(Country_pools[,i,,time+1])>0
          distance_array[,i,,]<-1/global_spread[3,g]
        }
        
        local_change<-log(colSums(Province_pools[,,,,time+1])/colSums(Province_pools[,,,,1]))
        country_change<-log(colSums(Country_pools[,,,time+1])/colSums(Country_pools[,,,1]))
        continent_change<-log(colSums(Regional_pools[,,time+1])/colSums(Regional_pools[,,1]))
        global_change<-log(sum(rowSums(Regional_pools[,,time+1])>0)/sum(rowSums(Regional_pools[,,1])>0))
        
        local_beta<-jaccard(betapart2(Province_pools[,,,,time+1],Province_pools[,,,,1]))
        country_beta<-jaccard(betapart2(Country_pools[,,,time+1],Country_pools[,,,1]))
        continent_beta<-jaccard(betapart2(Regional_pools[,,time+1], Regional_pools[,,1]))
        global_beta<-jaccard(betapart(rowSums(Regional_pools[,,time+1])>0, rowSums(Regional_pools[,,1])>0))
        
        Change.df<-rbind(Change.df,data.frame(Prob_local_loss=overall_extirp_prob,
                                              Prob_invasion=invasion_prob,
                                              Global_spread=global_spread[3,g],
                                              Scale_size=c(rep(1,length(local_change)),rep(Provinces_per_country,length(country_change)),rep(Provinces_per_country*Countries_per_region, length(continent_change)),Provinces_per_country*Countries_per_region*Regions),
                                              Scale=factor(c(rep("Local",length(local_change)),rep("Country",length(country_change)),rep("Continent", length(continent_change)),"Global"),levels=c("Local","Country","Continent","Global"),ordered = T),
                                              log_ratio=c(local_change,country_change,continent_change,global_change),
                                              jaccard=c(local_beta,country_beta,continent_beta,global_beta),
                                              time=time))
      }
    }
    Species.df<- rbind(Species.df, data.frame(richness = c(mean(colSums(Province_pools[,,,,1])),
                                                           mean(colSums(Country_pools[,,,1])), 
                                                           mean(colSums(Regional_pools[,,1])), 
                                                           mean(sum(rowSums(Regional_pools[,,1])>0)),
                                                           mean(colSums(Province_pools[,,,,time+1])),
                                                           mean(colSums(Country_pools[,,,time+1])), 
                                                           mean(colSums(Regional_pools[,,time+1])), 
                                                           mean(sum(rowSums(Regional_pools[,,time+1])>0))
    ),
    time = c(rep(1,4), rep(100, 4)),
    scale = c(1,10,100,1000),
    Prob_local_loss=overall_extirp_prob,
    Prob_invasion=invasion_prob,
    Global_spread=global_spread[3,g]))
  }
}

Change.df<-Change.df%>%
  mutate(Global_spread_text=Global_spread)%>%
  mutate(Global_spread_text=replace(Global_spread_text,Global_spread_text==1, "Global dispersal "),
         Global_spread_text=replace(Global_spread_text,Global_spread_text==100, "Intermediate dispersal"),
         Global_spread_text=replace(Global_spread_text,Global_spread_text==1000, "Localized dispersal"))%>%
  mutate(Prob_local_loss_text=paste("Loss =",Prob_local_loss))%>%
  mutate(Prob_invasion_text=paste("Colonization =",Prob_invasion))%>%
  mutate(Local_loss = Prob_local_loss_text) %>% 
  mutate(Local_loss = replace(Local_loss, Prob_local_loss == 0.01, "rare extirpation"),
         Local_loss = replace(Local_loss, Prob_local_loss == 0.03, "intermediate extirpation"),
         Local_loss = replace(Local_loss, Prob_local_loss == 0.1, "frequent extirpation")) %>% 
  mutate(Local_loss = factor(Local_loss, levels = rev(c("frequent extirpation", "intermediate extirpation", "rare extirpation")),ordered = T)) %>% 
  as.data.frame()

save(Change.df,file = "./workspace/non_equil.RData")
