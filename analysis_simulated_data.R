rm(list=ls())
library(magic)
library(igraph)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(vegan)
source('~/Dropbox/EAWAG PostDoc/06_scale_of_EWS/dispKernels.R',echo=F)
source('~/Dropbox/EAWAG PostDoc/06_scale_of_EWS/buildlandscape.R', echo=F)
source('~/Dropbox/EAWAG PostDoc/06_scale_of_EWS/Functions.R', echo=F)
source('~/Dropbox/Zurich PhD Research/1_Chapter_1/rolling_GAMs_methods.R', echo=F)

library(NetIndices)
library(betalink)

library(igraph)
library(DescTools)



nestedness_NODF <- function(web){
  web[web > 0] = 1
  SA <- nrow(web)
  SP <- ncol(web)
  N <- t(web) %*% web
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SP,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SP)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele == 0] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)] = 1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n1 <- sum(nes)
  
  N <- web %*% t(web)
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SA,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SA)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele ==0 ] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)]=1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n2 <- sum(nes)
  out <- 2*(n1 + n2) / (SA*(SA-1)+SP*(SP-1))
  return(out)
}



collapse.point<-function(dat,mut_strength,type_of_collapse, patches,species,dispersal_rate,
                         seed, webs){
  
  
  temp<-dat %>% filter(type.of.collapse == "global",
                       no_of_patches==patches,
                       no.of.species == species,
                       a==dispersal_rate,
                       random_seed==seed)
  
  rate.of.change.strength<- 0.15
  db<-numeric()
  for(r in 1:length( (temp$metacomm_biomass) -1)){
    
    db[r] <- (temp$metacomm_biomass[r+1] - temp$metacomm_biomass[r])
  }
  db[is.na(db)]<-0
  
  #450*rate.of.change.strength
  
  tipping.point.index<- max(which( abs(db) > 40 ))
  if(tipping.point.index == Inf) { tipping.point.index <- temp$g0[max(which(temp$metacomm_biomass <1 )) ]}
  else {tipping.point <- temp$g0[tipping.point.index]}
  
  
  if( max(abs(db)) > 40) {  collapse =1
  }else if(max(abs(db))< 40) {collapse = 0
  }
  
  
  point.of.collapse<- tipping.point
  #if(point.of.collapse == -Inf ){ 
   # point.of.collapse <- 0}
  
  
  return(list(point.of.collapse=point.of.collapse,
              #nestedness =temp$nestedness[1],
              network.size = species,
              #connectance=temp$connectance[1],
              collapse = collapse
  ))
  
}
aucsignal<-function(dat,mut_strength,type_of_collapse, patches, species ,dispersal_rate,
                    seed, tipping.point){
  
  temp<-dat %>% filter(type.of.collapse == "global",
                       no_of_patches== patches, 
                       no.of.species==species, 
                       a==dispersal_rate,
                       random_seed==seed)
  
  auc_ar1<-auc_sd<-auc_spatial_sd<-auc_spatial_variation<-auc_spatial_ar1<-auc_alpha_variation<-auc_spatial_correlation<-auc_regional_variability<-numeric()
  
  
  #temp<-te %>% filter(random_seed ==unique(te$random_seed)[i] )
  # temp<-na.omit(temp)
  if(tipping.point == 0){
    index_tp<-length(temp$g0)
  }else { index_tp <- length(temp$g0)
  }
  #area under the curve for species level ar1
  auc_ar1<-AUC(temp$g0[1:index_tp],(temp$species.ar1[1:index_tp]))
  
  #area under the curve for spatial autocorrelation ar1
  auc_spatial_ar1<-AUC(temp$g0[1:index_tp],(temp$spatial.ar1[1:index_tp]))
  
  
  #area under the curve for species level sd
  auc_sd<-AUC(temp$g0[1:index_tp], temp$species.sd[1:index_tp])
  
  #area under the curve for spatial  sd
  auc_spatial_sd<-AUC(temp$g0[1:index_tp], temp$spatial.sd[1:index_tp])
  
  #area under the curve for spatial variability
  auc_beta_variability<-AUC(temp$g0[1:index_tp], temp$beta_variability[1:index_tp])
  
  #area under the curve for alpha variaibility
  auc_alpha_variation<-AUC(temp$g0[1:index_tp], temp$alpha_variability[1:index_tp])
  
  #area under the curve for local level synchrony
  #auc_alpha_variation<-AUC(temp$g0[1:index_tp], temp$local_synchrony[1:index_tp])
  
  #area under the curve for spatial correlation
  auc_spatial_correlation<-AUC(temp$g0[1:index_tp], temp$spatial.correlation[1:index_tp])
  
  #area under the curve for regional variability
  auc_regional_variability<-as.numeric(AUC(temp$g0[1:index_tp], temp$metacommunity_variability[1:index_tp]))
  
  auc_local_synchrony <-as.numeric(AUC(temp$g0[1:index_tp], temp$local_synchrony[1:index_tp]))
  
  auc_regional_synchrony <- as.numeric(AUC(temp$g0[1:index_tp], temp$regional_synchrony[1:index_tp]))
  
  return(list(auc_regional_variability=auc_regional_variability,
              auc_spatial_correlation=auc_spatial_correlation,
              auc_alpha_variation=auc_alpha_variation,
              auc_beta_variability=auc_beta_variability,
              auc_spatial_sd=auc_spatial_sd,
              auc_spatial_ar1=auc_spatial_ar1,
              auc_sd=auc_sd,
              auc_ar1=auc_ar1,
              auc_local_synchrony=auc_local_synchrony,
              auc_regional_synchrony=auc_regional_synchrony,
              nestedness=temp$nestedness[1],
              network.size = temp$network.size[1],
              connectance = temp$connectance[1]
  ))
  
}



load(file = "Scale_of_EWS.RData")


#point of collapse estimation
point_of_collapse_data<-expand.grid(`a`= c(0,0.05,0.15),
                                    `no.of.species`=c(2,6,10),
                                    `no_of_patches` =c(2,5,15,20),
                                    `type.of.collapse` = "global",
                                    `random_seed`=4327+(1:1)*100) %>% 
  as_tibble %>% 
  mutate(`point_of_collaspe`=0,
         `collapse` = 0)

for(i in 1:nrow(point_of_collapse_data)){
  tempp<-collapse.point(dat = fact, mut_strength = fact$g0, 
                        type_of_collapse = point_of_collapse_data$type.of.collapse[i], 
                        patches = point_of_collapse_data$no_of_patches[i],
                        species = point_of_collapse_data$no.of.species[i], 
                        dispersal_rate = point_of_collapse_data$a[i] , 
                        seed = point_of_collapse_data$random_seed[i])
  
  point_of_collapse_data$point_of_collaspe[i]<-tempp$point.of.collapse
  point_of_collapse_data$collapse[i]<-tempp$collapse
}



load("Species_tipping_points_spatial_data.RData")
#plotting point of collapse
tempdat<-point_of_collapse_data #%>% filter(no_of_patches !=20)

tempdat$point_of_collaspe<-tempdat$point_of_collaspe
hist(tempdat$point_of_collaspe)
qqnorm((tempdat$point_of_collaspe)^4)
qqline(tempdat$point_of_collaspe,
       col="red")

summary(model1<-glm(point_of_collaspe/5 ~ no.of.species*no_of_patches+a,
                    family = quasibinomial,
                    data=tempdat))

#rule of thumb of over dispersion: the ration of residual deviance to df should be 1. 


point_of_collapse_data %>%  
  filter_all(all_vars(!is.infinite(.))) %>% 
  # filter(no_of_patches !=20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = no.of.species,
             y=point_of_collaspe, 
             color = factor(no_of_patches)))+
  #geom_point(outlier.size = 0)+
  geom_boxplot(alpha = 5/10, size = 1)+
  ylab("Point of transition")+
  xlab("Network size")+
 # stat_smooth(method = "glm",se = T,
  #            method.args = list(family = "quasibinomial"))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)



#area under the curve for signals of collapse

auc_ews_data<-expand.grid(`a`= c(0,0.05,0.15),
                          `no.of.species`=c(2,6,10),
                          `no_of_patches` =c(2,5,10,20),
                          `type.of.collapse` =  "global",
                          `random_seed`=4327+(1:1)*100) %>% 
  as_tibble %>% 
  mutate(`auc_regional_variability`=0,
         `auc_spatial_correlation`=0,
         `auc_alpha_variation`=0,
         `auc_spatial_variation`=0,
         `auc_spatial_sd`=0,
         `auc_sd`=0,
         `auc_ar1`=0,
         `spatial_ar1`=0)


for(i in 1:nrow(auc_ews_data)){
  dtt<-aucsignal(dat = fact, 
                 mut_strength = fact$g0,
                 type_of_collapse =auc_ews_data$type.of.collapse[i],
                 patches = auc_ews_data$no_of_patches[i],
                 species = auc_ews_data$no.of.species[i],
                 dispersal_rate = auc_ews_data$a[i],
                 seed = auc_ews_data$random_seed[i],
                 tipping.point = point_of_collapse_data$point_of_collaspe[i])
  
  
  auc_ews_data$auc_regional_variability[i] <-dtt$auc_regional_variability
  auc_ews_data$auc_spatial_correlation[i] <- dtt$auc_spatial_correlation
  auc_ews_data$auc_spatial_variation[i]<-dtt$auc_spatial_variation
  auc_ews_data$auc_spatial_sd[i] <-dtt$auc_spatial_sd
  auc_ews_data$auc_sd[i] <-dtt$auc_sd
  auc_ews_data$auc_ar1[i] <- dtt$auc_ar1
  auc_ews_data$auc_alpha_variation[i] <-dtt$auc_alpha_variation
  auc_ews_data$spatial_ar1[i] <- dtt$auc_spatial_ar1
  
}

# auc ews 
tempauc<-auc_ews_data #%>% filter(no_of_patches !=20)

summary(model1<-lm(log(auc_spatial_variation) ~ no.of.species*no_of_patches*a,
                   data=tempauc))
hist(log(tempauc$auc_spatial_variation))
plot(model1)

auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  #filter(no_of_patches != 20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = no.of.species,
             y=auc_spatial_variation, 
             color = factor(no_of_patches)))+
  #geom_boxplot(outlier.size = 0)+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC spatial deviation ")+
  xlab("Network size")+
  stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)



# auc standard deviation
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  #filter(no_of_patches != 20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = no.of.species,
             y=auc_sd, 
             color = factor(no_of_patches)))+
  #geom_boxplot(outlier.size = 0)+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC species level standard deviation ")+
  xlab("Network size")+
  # stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)


#auc alpha variability
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  #filter(no_of_patches != 20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = no.of.species,
             y=auc_alpha_variation, 
             color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC alpha variation")+
  xlab("Network size")+
  #stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)
#auc_regional_variability


auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  filter(type.of.collapse ==  "global") %>% 
  #filter(no_of_patches != 20) %>% 
  ggplot(aes(x = no.of.species,
             y=auc_regional_variability, 
             color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC regional variation")+
  xlab("Network size")+
  stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)

#spatial autocorrelation at first lag1
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  ggplot(aes(x = no.of.species,
             y=spatial_ar1, 
             color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC regional variation")+
  xlab("Network size")+
  #stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)


#spatial autocorrelation at first lag1
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  ggplot(aes(x = no.of.species,
             y=auc_ar1, 
             color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC regional variation")+
  xlab("Network size")+
  #stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)

