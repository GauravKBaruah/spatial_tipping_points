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



load("scale_of_ews_webs_all_data.RData")
#fact$nestedness<-round(fact$nestedness,4)

mydir = 'plant_pollinator'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
myfiles<-myfiles[1:56]



collapse.point<-function(dat,mut_strength,type_of_collapse, patches ,dispersal_rate,
                         seed, webs){
  
  
  temp<-dat %>% filter(type.of.collapse == type_of_collapse,
                       no_of_patches==patches,
                       web == webs,
                       a==dispersal_rate,
                       random_seed==seed)
  
  rate.of.change.strength<- 0.15
  db<-numeric()
  for(r in 1:length( (temp$metacomm_biomass) -1)){
    
    db[r] <- (temp$metacomm_biomass[r+1] - temp$metacomm_biomass[r])
  }
  db[is.na(db)]<-0
  
  #450*rate.of.change.strength
  
  tipping.point.index<-  min(which(temp$metacomm_biomass < (temp$metacomm_biomass[1] -0.95*temp$metacomm_biomass[1])))  #max(which( abs(db) > 50 ))
  tipping.point <- temp$g0[tipping.point.index]
 
 
  if( max(abs(db)) > 20) {  collapse =1
  }else if(max(abs(db))< 20) {collapse = 0
  }
  
  
  point.of.collapse<- tipping.point
 # if(point.of.collapse == Inf ){ 
  #  point.of.collapse <- NA }
#  else (point.of.collapse <-tipping.point )
  
  
  return(list(point.of.collapse=point.of.collapse,
              nestedness =temp$nestedness[1],
              network.size = temp$network.size[1],
              connectance=temp$connectance[1],
              collapse = collapse
              ))
  
}
aucsignal<-function(dat,mut_strength,type_of_collapse, patches, species ,dispersal_rate,
                    seed, tipping.point){
  
  temp<-dat %>% filter(type.of.collapse == type_of_collapse,
                       no_of_patches== patches, 
                       network.size==species, 
                       a==dispersal_rate,
                       random_seed==seed)
  
  auc_ar1<-auc_sd<-auc_spatial_sd<-auc_spatial_variation<-auc_spatial_ar1<-auc_alpha_variation<-auc_spatial_correlation<-auc_regional_variability<-numeric()
  
  
  #temp<-te %>% filter(random_seed ==unique(te$random_seed)[i] )

  index_tp<-length(temp$g0)
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




#point of collapse estimation
point_of_collapse_data<-expand.grid(`a`= c(0,0.05,0.15),
                                    `web`=myfiles,
                                    #`nestedness`=unique(fact$nestedness),
                                    `no_of_patches` =c(2,5,10,20),
                                    `type.of.collapse` = "global",
                                    `random_seed`=4327+(1:1)*100) %>% 
  as_tibble %>% 
  mutate(`point_of_collaspe`=0,
         `network.size`=0,
         `connectance`=0,
         `nestedness`=0,
         `collapse`=0)

for(i in 1:nrow(point_of_collapse_data)){
  tempp<-collapse.point(dat = fact, 
                        mut_strength = fact$g0, 
                        type_of_collapse = point_of_collapse_data$type.of.collapse[i], 
                        patches = point_of_collapse_data$no_of_patches[i],
                        webs = as.character(point_of_collapse_data$web[i][1]), 
                        dispersal_rate = point_of_collapse_data$a[i] , 
                        seed = point_of_collapse_data$random_seed[i])
  
  point_of_collapse_data$network.size[i] <-tempp$network.size
  point_of_collapse_data$nestedness[i] <-tempp$nestedness
  point_of_collapse_data$connectance[i] <-tempp$connectance
  point_of_collapse_data$point_of_collaspe[i]<-tempp$point.of.collapse
  point_of_collapse_data$collapse[i]<- tempp$collapse

}




#plotting point of collapse
tempdat<-point_of_collapse_data #%>% filter(no_of_patches !=20)

tempdat$point_of_collaspe<-tempdat$point_of_collaspe
hist((log(tempdat$point_of_collaspe)))
qqnorm((tempdat$point_of_collaspe))

summary(model1<-glm((point_of_collaspe)/5 ~ network.size*no_of_patches+nestedness + connectance,
                    family="quasibinomial",
                   data=tempdat))

summary(model2<-glm( (point_of_collaspe)/5 ~ nestedness*no_of_patches+a,
                     family="quasibinomial",
                     data=tempdat))

summary(model2<-lm(log(point_of_collaspe) ~ connectance*no_of_patches+a,
                    #family = quasibinomial,
                    data=tempdat))


#rule of thumb of over dispersion: the ration of residual deviance to df should be 1. 

formula <- y ~ poly(x, 2, raw = TRUE)

point_of_collapse_data %>%  
 # filter_all(all_vars(!is.infinite(.))) %>% 
  # filter(no_of_patches !=20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = network.size,
             y=point_of_collaspe, 
             color = factor(no_of_patches)))+
  #geom_point(outlier.size = 0)+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 1,outlier.size = 0)+
  ylab("Point of transition")+
  xlab("Network size")+
  stat_smooth(method = "lm",se = T,   method.args = list(family = "normal"))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)


point_of_collapse_data %>%  
  # filter_all(all_vars(!is.infinite(.))) %>% 
  # filter(no_of_patches !=20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = nestedness,
             y=point_of_collaspe/5, 
             color = factor(no_of_patches)))+
  #geom_point(outlier.size = 0)+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("Point of transition")+
  xlab("Nestedness")+
  stat_smooth(method = "glm",se = F,
              method.args = list(family = "quasibinomial"))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)


point_of_collapse_data %>%  
  # filter_all(all_vars(!is.infinite(.))) %>% 
  # filter(no_of_patches !=20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = connectance,
             y=point_of_collaspe, 
             color = factor(no_of_patches)))+
  #geom_point(outlier.size = 0)+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("Point of transition")+
  xlab("Connectance")+
  stat_smooth(method = "lm",se = F,
              method.args = list(family = "normal"))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)

#area under the curve for signals of collapse

auc_ews_data<-expand.grid(`a`= c(0,0.05,0.15),
                          `no.of.species`=as.numeric(levels((as.factor(fact$network.size)))),
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
         `spatial_ar1`=0,
         `nestedness` = 0,
         `network.size`=0,
         `connectance`=0,
         auc_beta_variability=0,
         auc_regional_synchrony=0,
         auc_local_synchrony=0
         )


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
  auc_ews_data$auc_beta_variability[i]<-dtt$auc_beta_variability
  auc_ews_data$auc_spatial_sd[i] <-dtt$auc_spatial_sd
  auc_ews_data$auc_sd[i] <-dtt$auc_sd
  auc_ews_data$auc_ar1[i] <- dtt$auc_ar1
  auc_ews_data$auc_alpha_variation[i] <-dtt$auc_alpha_variation
  auc_ews_data$spatial_ar1[i] <- dtt$auc_spatial_ar1
  auc_ews_data$auc_regional_synchrony[i] <- dtt$auc_regional_synchrony
  auc_ews_data$auc_local_synchrony[i] <- dtt$auc_local_synchrony
  auc_ews_data$nestedness[i] <- dtt$nestedness
  auc_ews_data$connectance[i]<- dtt$connectance
  auc_ews_data$network.size[i] <-dtt$network.size
  
}

# auc ews 
tempauc<-auc_ews_data #%>% filter(no_of_patches !=20)

summary(model1<-lm(log(auc_spatial_variation) ~ no.of.species*no_of_patches*a,
                   data=tempauc))
hist(log(tempauc$auc_spatial_variation))
plot(model1)


#beta variability

auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  #filter(no_of_patches != 20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = network.size,
             y=auc_beta_variability, 
             color = factor(no_of_patches)))+
  #geom_boxplot(outlier.size = 0)+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC spatial variability ")+
  xlab("Network size")+
  stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)



# auc patch variability
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  #filter(no_of_patches != 20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = no.of.species,
             y=auc_alpha_variation, 
             color = factor(no_of_patches)))+
  #geom_boxplot(outlier.size = 0)+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21, outlier.size =0)+
  ylab("AUC alpha variability ")+
  xlab("Network size")+
   stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)


#auc metacommunity variability
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  #filter(no_of_patches != 20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = network.size,
             y=auc_regional_variability, 
             color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC regional variability")+
  xlab("Network size")+
 stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)



#auc metacommunity syncrhony
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  #filter(no_of_patches != 20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = network.size,
             y=auc_regional_synchrony, 
             color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC metacommunity syncrhony")+
  xlab("Network size")+
  stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)


#auc local syncrhony
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  #filter(no_of_patches != 20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = network.size,
             y=auc_local_synchrony, 
             color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC local syncrhony")+
  xlab("Network size")+
  stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)


#spatial SD
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  ggplot(aes(x = no.of.species,
             y=auc_spatial_sd, 
             color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC spatial variation")+
  xlab("Network size")+
  stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)


# auc SD at species level

auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  ggplot(aes(x = no.of.species,
             y=auc_sd, 
             color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC regional variation")+
  xlab("Network size")+
  stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)


#auc species level ar1
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
ggplot(aes(x = no.of.species,
           y=auc_ar1, 
           color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC regional variation")+
  xlab("Network size")+
 # stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)




#spatial autocorrelation at first lag1
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  ggplot(aes(x = network.size,
             y=auc_ar1, 
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
  ggplot(aes(x = network.size,
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
  ggplot(aes(x = network.size,
             y=auc_spatial_correlation, 
             color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC regional variation")+
  xlab("Network size")+
  #stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)






fact$network.size_class <- findInterval(fact$network.size, seq(10, max(fact$network.size), 12), left.open = TRUE)

#biomass
fact %>% 
  ggplot(aes(x=g0, y=alpha.richness, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("metacommunity biomass")+xlab("mutualistic strength")+
  geom_smooth(method = 'loess',se=TRUE)+theme_bw()+
  facet_grid(a~network.size_class)


#species level ar1
fact %>%  filter(type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=species.ar1, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("species ar1")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~network.size_class)


fact %>%  filter(type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=alpha_variability, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("patch variability")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~network.size_class)

# spatial Standard deviation

fact %>%  filter( type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=spatial.sd, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("spatial SD")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~network.size_class)




#spatial correlation


fact %>%  filter(type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=spatial.correlation, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("spatial correlation")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~network.size_class)




# regional synchrony


fact %>%  filter( type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=regional_synchrony, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("Regional synchrony")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~network.size_class)



fact %>%  filter( type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=local_synchrony, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("local synchrony")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~network.size_class)



