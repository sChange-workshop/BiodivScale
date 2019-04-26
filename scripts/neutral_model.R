### HEADER #####################################################################
##' @title Equilibirum neutral biodiversity change simulation
##'
##' @author Patrick Thompson
##' @contact patrick.thompson@zoology.ubc.ca
##' 
##' @date 2017
##' 
##' @description Simulates scale dependent biodiversity change from local to global scales. 
##' Based on Hubbell's neutral model
##' 
##' 
##' ################################################################################

library(dplyr)
library(tidyr)

neutral.metacom_change<-function(Dv = c(5,10), MC_dispersalV = c(0.1,0.1), gens = 1000, abundanceV = c(50, 50), MC_size = 10, scale_3_size = 10, scale_4_size = 3, spec_prob = 0.005, initialize = 1){
  D<-Dv[1]
  MC_dispersal<-MC_dispersalV[1]
  abundance<-abundanceV[1]
  change<-FALSE
  
  if(initialize == 1){
    comA<-data.frame(scale4 = rep(1:scale_4_size,each = MC_size*abundance*scale_3_size),scale3 = rep(1:scale_3_size,each = MC_size*abundance), site = rep(1:MC_size,each = abundance),species = rep(1:(MC_size*scale_3_size*scale_4_size),each = abundance))
  } else {
    comA<-data.frame(scale4 = rep(1:scale_4_size,each = MC_size*abundance*scale_3_size),scale3 = rep(1:scale_3_size,each = MC_size*abundance), site = rep(1:MC_size,each = abundance),species = 1)
  }
  diversity.df<-data.frame()
  pb <- txtProgressBar(min = 0, max = gens*2, style = 3)
  for(i in 1:gens){
    maxSp<-max(comA$species)
    comA<-comA %>% 
      group_by(site, scale3, scale4) %>% 
      mutate(survive = sample(c(rep(FALSE,D),rep(TRUE, abundance - D)),size = abundance,replace = FALSE)) %>%
      mutate(species = replace(species, survive == FALSE, sample(species[survive == TRUE],size = D, replace = TRUE))) %>%
      mutate(disperse = FALSE) %>%
      mutate(disperse = replace(disperse, survive == FALSE, rbinom(n = D,size = 1,prob = MC_dispersal))) %>% 
      ungroup() %>% 
      group_by(scale4) %>% 
      mutate(disperse3 = FALSE) %>% 
      mutate(disperse3 = replace(disperse3, disperse == TRUE, rbinom(n = sum(disperse),size = 1,prob = MC_dispersal))) %>% 
      ungroup() %>% 
      mutate(disperse4 = FALSE) %>% 
      mutate(disperse4 = replace(disperse4, disperse3 == TRUE, rbinom(n = sum(disperse3),size = 1,prob = MC_dispersal))) %>% 
      group_by(site) %>% 
      mutate(species = replace(species, disperse == TRUE, sample(species,size = sum(disperse), replace = TRUE))) %>% 
      ungroup() %>% 
      mutate(species = replace(species, disperse3 == TRUE, sample(species,size = sum(disperse3), replace = TRUE))) %>% 
      mutate(speciate = FALSE) %>% 
      mutate(speciate = replace(speciate,survive == F, rbinom(n = MC_size*D,size = 1,prob = spec_prob))) %>% 
      mutate(species = replace (species, speciate == TRUE, c(1:sum(speciate)+maxSp)))
    
    diversity.df<-bind_rows(diversity.df,bind_rows(comA %>% 
                                                     group_by(site,scale3,scale4) %>% 
                                                     summarise(S = length(unique(species))) %>% 
                                                     mutate(scale = "1") %>% 
                                                     mutate(id = paste(1,scale4, scale3,site, sep = "_")) %>% 
                                                     ungroup() %>% 
                                                     select(S,scale,id),
                                                   comA %>% 
                                                     group_by(scale3,scale4) %>% 
                                                     summarise(S = length(unique(species))) %>% 
                                                     mutate(scale = "2") %>% 
                                                     mutate(id = paste(1,scale4, scale3, sep = "_")) %>% 
                                                     ungroup() %>% 
                                                     select(S,scale,id),
                                                   comA %>% 
                                                     group_by(scale4) %>% 
                                                     summarise(S = length(unique(species))) %>% 
                                                     mutate(scale = "3") %>% 
                                                     mutate(id = paste(1,scale4, sep = "_")) %>% 
                                                     ungroup() %>% 
                                                     select(S,scale,id),
                                                   comA %>% 
                                                     summarise(S = length(unique(species))) %>% 
                                                     mutate(scale = "4") %>% 
                                                     mutate(id = paste(1, sep = "_")) %>% 
                                                     select(S,scale,id)) %>% 
                              mutate(time = i, change = FALSE))
    setTxtProgressBar(pb, i)
  }
  
  D<-Dv[2]
  MC_dispersal<-MC_dispersalV[2]
  abundance<-abundanceV[2]
  
  comA<-comA %>% 
    group_by(site,scale3,scale4) %>%
    filter(row_number()<=abundance)
  
  if(abundanceV[2]>abundanceV[1]){
    difference<-abundanceV[2]-abundanceV[1]
    comA<-bind_rows(comA,comA %>% 
                      group_by(site, scale3, scale4) %>% 
                      filter(row_number()<=10))
  }
  
  for(i in (gens+1):(gens*2)){
    maxSp<-max(comA$species)
    comA<-comA %>% 
      group_by(site, scale3, scale4) %>% 
      mutate(survive = sample(c(rep(FALSE,D),rep(TRUE, abundance - D)),size = abundance,replace = FALSE)) %>%
      mutate(species = replace(species, survive == FALSE, sample(species[survive == TRUE],size = D, replace = TRUE))) %>%
      mutate(disperse = FALSE) %>%
      mutate(disperse = replace(disperse, survive == FALSE, rbinom(n = D,size = 1,prob = MC_dispersal))) %>% 
      ungroup() %>% 
      group_by(scale4) %>% 
      mutate(disperse3 = FALSE) %>% 
      mutate(disperse3 = replace(disperse3, disperse == TRUE, rbinom(n = sum(disperse),size = 1,prob = MC_dispersal))) %>% 
      ungroup() %>% 
      mutate(disperse4 = FALSE) %>% 
      mutate(disperse4 = replace(disperse4, disperse3 == TRUE, rbinom(n = sum(disperse3),size = 1,prob = MC_dispersal))) %>% 
      group_by(site) %>% 
      mutate(species = replace(species, disperse == TRUE, sample(species,size = sum(disperse), replace = TRUE))) %>% 
      ungroup() %>% 
      mutate(species = replace(species, disperse3 == TRUE, sample(species,size = sum(disperse3), replace = TRUE))) %>% 
      mutate(speciate = FALSE) %>% 
      mutate(speciate = replace(speciate,survive == F, rbinom(n = MC_size*D,size = 1,prob = spec_prob))) %>% 
      mutate(species = replace (species, speciate == TRUE, c(1:sum(speciate)+maxSp)))
    
    diversity.df<-bind_rows(diversity.df,bind_rows(comA %>% 
                                                     group_by(site,scale3,scale4) %>% 
                                                     summarise(S = length(unique(species))) %>% 
                                                     mutate(scale = "1") %>% 
                                                     mutate(id = paste(1,scale4, scale3,site, sep = "_")) %>% 
                                                     ungroup() %>% 
                                                     select(S,scale,id),
                                                   comA %>% 
                                                     group_by(scale3,scale4) %>% 
                                                     summarise(S = length(unique(species))) %>% 
                                                     mutate(scale = "2") %>% 
                                                     mutate(id = paste(1,scale4, scale3, sep = "_")) %>% 
                                                     ungroup() %>% 
                                                     select(S,scale,id),
                                                   comA %>% 
                                                     group_by(scale4) %>% 
                                                     summarise(S = length(unique(species))) %>% 
                                                     mutate(scale = "3") %>% 
                                                     mutate(id = paste(1,scale4, sep = "_")) %>% 
                                                     ungroup() %>% 
                                                     select(S,scale,id),
                                                   comA %>% 
                                                     summarise(S = length(unique(species))) %>% 
                                                     mutate(scale = "4") %>% 
                                                     mutate(id = paste(1, sep = "_")) %>% 
                                                     select(S,scale,id)) %>% 
                              mutate(time = i, change = TRUE))
    setTxtProgressBar(pb, i)
  }
  
  
  
  
  
  close(pb)
  return(diversity.df)
}

A <- c(5)
B <- c(0.1, 0.2, 0.05)
D <- c(50, 40, 60)

LRR<-data.frame()
for(reps in 1:10){
  for(a in A){
    for(b in B){
      for(d in D){
        model<-neutral.metacom_change(Dv = c(5,a),MC_dispersalV = c(0.1,b),abundanceV = c(50,d),gens = 1000,initialize = 1)
        
        #print(ggplot(model, aes(x= time, y = S, color = scale, group = id))+
        #  geom_line()+
        #  scale_y_log10())
        
        LRR<-bind_rows(LRR,model %>% 
                         #filter(time %in% seq(1000,2000, by=100)) %>% 
                         group_by(scale,id) %>%
                         summarise(LRR = log(mean(S[time <= 2000 & time>1900])/mean(S[time <= 1000 & time>900]))) %>% 
                         mutate(death_change = 5-a, dispersal_change = 0.1-b, abundance_change = 50-d, rep = reps))
      }
    }
  }  
}

save(LRR, file = "./workspace/neutral_model_out.RData")