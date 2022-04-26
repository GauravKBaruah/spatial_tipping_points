rm(list=ls())
library(magic)
library(igraph)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(vegan)
source('~/Dropbox/EAWAG PostDoc/06_scale_of_EWS/spatial_tipping_points/dispKernels.R',echo=F)
source('~/Dropbox/EAWAG PostDoc/06_scale_of_EWS/spatial_tipping_points/buildlandscape.R', echo=F)
source('~/Dropbox/EAWAG PostDoc/06_scale_of_EWS/spatial_tipping_points/Functions.R', echo=F)


library(NetIndices)
library(betalink)
library(igraph)
library(DescTools)

create_matrix<-function(SA,SP,required_nestedness){

nestd<-0

while(nestd !=required_nestedness) {
  
  web <- matrix(rbinom(SA*SP, 1, prob=0.6),nrow=SA,ncol =SP)
  nestd<-round(nestedness_NODF(web),1)
  
  if(nestd == required_nestedness){
    web<-web
  }else {
    nestd <- nestd
  }
  
}
return(web)
  
}


#empty arrays for abundance of species
start.time <-2000


#coordinates for the 2-D patches
d<-2 #no. of dimensions
L<-2 #length
coordinates_2  <- matrix(runif(2 * d), 2, d) * L # random coordinates for N patches in d dimensions
coordinates_5  <- matrix(runif(5 * d), 5, d) * L # random coordinates for N patches in d dimensions
coordinates_10 <- matrix(runif(10 * d), 10, d) * L # random coordinates for N patches in d dimensions
coordinates_20 <- matrix(runif(20 * d), 20, d) * L # random coordinates for N patches in d dimensions
Cordinates_list<-list(coordinates_2=coordinates_2,
                      coordinates_5=coordinates_5,
                      coordinates_10=coordinates_10,
                      coordinates_20=coordinates_20)



patch_removed_2<-sample(1:2,1,replace = F)
patch_removed_5<-sample(1:5,1,replace = F)
patch_removed_10<-sample(1:10,1,replace = F)
patch_removed_20<-sample(1:20,1,replace = F)


Disp.mat<- BuildRndLandscape(N = 2,w=0.5, Coordinates=coordinates_2,FUN = thompson_etal, sigmaA = 0,Bounds = 'open')
plot_net(Disp.mat,coordinates_2)


model<-"Thompson_et_al"
pars<-list()


mydir = 'plant_pollinator'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
myfiles<-myfiles[1:56]


#creating an empty data frame
fact<- expand.grid(`w`=1,
                   `g0` = seq(5,0,-0.175),
                   `a`= c(0,0.05,0.15),
                   `web`=myfiles,
                   `no_of_patches` =c(2,5,10,20),
                   `type.of.collapse` =  "global",
                   `random_seed`=4327+(1:1)*100) %>%
  as_tibble %>%
  mutate(`species.ar1`=0,
         `spatial.ar1`=0,
         `spatial.sd`=0,
         `beta_variability`=0,
         `alpha_variability`=0,
         `metacommunity_variability`=0,
         `metacomm_biomass`=0,
         `local_synchrony`=0,
         `regional_synchrony`=0,
         `network_dissimilarity`= 0,
         `spatial.correlation`=0,
         `species.ar1`=0,
          `species.sd` = 0,
         `patch_affected`=0,
         `nestedness`=0,
         `connectance`= 0,
         `network.size`=0,
         `mean.biomass.patch.affected`=0,
         `patch_richness`=0,
         `metacommunity_richness`=0,
         `beta_richness` =0)

for(i in 1:nrow(fact)){
  
  g<-adj.mat(myfiles[which(myfiles == fact$web[i])]) #network web names
  #g<-g[-1,-1] 
  fact$network.size[i]<- nrow(g)+ncol(g)
  fact$nestedness[i] <- nestedness_NODF(g)
  
}




for(j in 1:nrow(fact)){

  if(fact$no_of_patches[j] == 2){
    coordinates<- Cordinates_list$coordinates_2
  }else if(fact$no_of_patches[j] == 5){
    coordinates<- Cordinates_list$coordinates_5
  }else if(fact$no_of_patches[j] == 10){
    coordinates<- Cordinates_list$coordinates_10
  }else if(fact$no_of_patches[j] == 20){
    coordinates<- Cordinates_list$coordinates_20
  }
  
  
  if(fact$no_of_patches[j] == 2 & fact$type.of.collapse[j] == "local"){
    coordinates<- Cordinates_list$coordinates_2
    patch_affected_2<-patch_removed_2
  }else if(fact$no_of_patches[j] == 5 & fact$type.of.collapse[j] == "local"){
    coordinates<- Cordinates_list$coordinates_5
    patch_affected_5<- patch_removed_5
  }else if(fact$no_of_patches[j] == 10 & fact$type.of.collapse[j] == "local"){
    coordinates<- Cordinates_list$coordinates_10
    patch_affected_15 <- patch_removed_15
  }else if(fact$no_of_patches[j] == 20 & fact$type.of.collapse[j] == "local"){
    coordinates<- Cordinates_list$coordinates_20
    patch_affected_20<- patch_removed_20
  }
  
  
  Disp.mat<-BuildRndLandscape(N = fact$no_of_patches[j],w=1,Coordinates=coordinates,
                              FUN = thompson_etal, sigmaA = 0,Bounds = 'open')
  
  
  g<-adj.mat(myfiles[which(myfiles == fact$web[j])]) 
  #g<-g[-1,-1]
  
  Aspecies<- nrow(g) # no of animal species
  Pspecies<- ncol(g) # no of plant species
  degree.animals<-degree.plants<-numeric()
  #adj_matrix<-adj.matr(g=g)
  
  alpha_a <- array(NA,dim=c(fact$no_of_patches[j], Aspecies, Aspecies ))
  alpha_p <-array(NA,dim=c(fact$no_of_patches[j], Pspecies, Pspecies ))
  ra <-  matrix(NA, fact$no_of_patches[j], Aspecies)
  rp <-  matrix(NA, fact$no_of_patches[j], Pspecies)
  gamma<-array(NA,dim=c(fact$no_of_patches[j], Aspecies, Pspecies ))
  

     
      for(k in 1:fact$no_of_patches[j]){
        alpha_a[k,,] <- runif(Aspecies^2, 0.005,0.01)
        alpha_p[k,,] <-runif(Pspecies^2,0.005,0.01)
       diag(alpha_a[k,,])<-0.5;diag(alpha_p[k,,])<-0.5
       
        # }else{ alpha_a <-alpha_a;alpha_p<-alpha_p
        # }
        # 
        
        ra[k,] <- runif(Aspecies,-0.1,-0.01)
        rp[k,] <- runif(Pspecies,-0.1,-0.01)
       if(fact$type.of.collapse[j] == "local" & fact$no_of_patches[j] == 2){  
        gamma[patch_affected_2,,] <- rnorm(Aspecies*Pspecies, fact$g0[j], 0.005)
       }else if(fact$type.of.collapse[j] == "local" & fact$no_of_patches[j] == 5) {
         gamma[patch_affected_5,,] <- rnorm(Aspecies*Pspecies, fact$g0[j], 0.005)
       }else if(fact$type.of.collapse[j] == "local" & fact$no_of_patches[j] == 15) {
         gamma[patch_affected_15,,] <- rnorm(Aspecies*Pspecies, fact$g0[j], 0.005)
       }else if(fact$type.of.collapse[j] == "local" & fact$no_of_patches[j] == 20) {
         gamma[patch_affected_20,,] <- rnorm(Aspecies*Pspecies, fact$g0[j], 0.005)
       }else if(fact$type.of.collapse[j] == "global"){
         for(i in 1:Aspecies){
           for(r in 1:Pspecies){
             gamma[k,i,r] <- rnorm(1,fact$g0[j], 0.00)*g[i,r]
           }
         }
       }
      }
      
        
  Ga <- fact$a[j]
  Gp <- fact$a[j]
  
      
  pars<-list(Aspecies=Aspecies,Pspecies=Pspecies,patches=fact$no_of_patches[j], Disp.mat=Disp.mat,Ga=Ga,Gp=Gp,g0=fact$g0[j],
             alpha_a=alpha_a,alpha_p=alpha_p,ra=ra,rp=rp,gamma=gamma,model="Thompson_et_al",type_of_collapse= fact$type.of.collapse[j])
  
  start.time = 4000
  state<-0.5
  model.t<-lapply(1, Mcommunity_1,time=start.time,state=state,
                  pars=pars)
  #model.t[[1]]$N[800,1,]
  
  
  if(fact$type.of.collapse[j] == "global") {
    fact$patch_affected[j] <- 0
  }else if (fact$type.of.collapse[j] == "local" & fact$no_of_patches[j] == 2) 
  {fact$patch_affected[j] <- patch_removed_2
  }else if (fact$type.of.collapse[j] == "local" & fact$no_of_patches[j] == 5) {
    fact$patch_affected[j] <- patch_removed_5
  }else if (fact$type.of.collapse[j] == "local" & fact$no_of_patches[j] == 10) {
    fact$patch_affected[j] <- patch_removed_10
  }else if (fact$type.of.collapse[j] == "local" & fact$no_of_patches[j] == 20) {
    fact$patch_affected[j] <- patch_removed_20
  }

  temp_mat<-meta.stats(dat=model.t, patches=pars$patches,
                       species=(pars$Aspecies+pars$Pspecies),time=start.time, 
                       disp.mat=Disp.mat,
                       type_of_collapse = fact$type.of.collapse[j],
                       patch_affected =  fact$patch_affected[j])
  
  #variability in the regional scales of the metacommunity
  fact$metacommunity_variability[j] <-temp_mat$metacom_variability
  fact$alpha_variability[j] <- temp_mat$alpha.variability
  fact$beta_variability[j]<- temp_mat$beta.variability
  fact$metacomm_biomass[j] <- temp_mat$metacomm_biomass
  fact$local_synchrony[j] <- temp_mat$avg.local.synchrony
  fact$regional_synchrony[j] <- temp_mat$reg.synchrony
  fact$spatial.ar1[j] <-temp_mat$Acf_patch
  
  fact$patch_richness[j] <-temp_mat$alpha.richness
  fact$metacommunity_richness[j] <-temp_mat$gamma.richness
  fact$beta_richness[j] <- temp_mat$spatial.beta.richness
  
  fact$network.size[j] <- (Aspecies+Pspecies)
  fact$nestedness[j] <- nestedness_NODF(g)
  fact$connectance[j] <- Connectance(g)
 
   fact$species.ar1[j]<-mean(temp_mat$species_level_acf,na.rm=T)
  fact$spatial.sd[j]<-temp_mat$metacom_variability
  fact$spatial.correlation[j]<-temp_mat$spatial.correlation
  fact$species.sd[j]<-mean(temp_mat$species_level_sd,na.rm=T)
   if(fact$type.of.collapse[j] == "global"){ 
     fact$mean.biomass.patch.affected[j] <- 0
   }else if(fact$type.of.collapse[j] == "local"){
     fact$mean.biomass.patch.affected[j] <- temp_mat$mean.biomass.patch.affected
     }

  print(j)
}

# #save(fact, file="example_network_collapse.RData")
# 
# #15019
# #
# #save(fact, file="scale_of_ews_real_webs_20may.RData")
# #load("scale_of_ews_real_webs_20may.RData")
# #load("scale_of_ews_real_webs_15may.RData")
# #fact<-rbind(fact,fact_1)
# 
# #save(fact,file="scale_of_ews_webs_all_data.RData")
# 
# load("scale_of_ews_webs_all_data.RData")
# #load("example_network_collapse.RData")
# 
# #point of collapse estimation
# 
# 
# fact<-fact %>% filter(web != "plant_pollinator/M_PL_013.csv" & web!= "plant_pollinator/M_PL_046.csv")
# 
# point_of_collapse_data<-expand.grid(`a`= c(0,0.05,0.15),
#                                     `no.of.species`=as.numeric(levels(as.factor(fact$network.size))),
#                                     `no_of_patches` =c(2,5,10,20),
#                                     `type.of.collapse` = "global",
#                                     `random_seed`=4327+(1:1)*100) %>% 
#   as_tibble %>% 
#   mutate(`point_of_collaspe`=0)
# 
# for(i in 1:nrow(point_of_collapse_data)){
#   tempp<-collapse.point(dat = fact, mut_strength = fact$g0, 
#                         type_of_collapse = point_of_collapse_data$type.of.collapse[i], 
#                         patches = point_of_collapse_data$no_of_patches[i],
#                         species = point_of_collapse_data$no.of.species[i], 
#                         dispersal_rate = point_of_collapse_data$a[i] , 
#                         seed = point_of_collapse_data$random_seed[i])
#   
#   point_of_collapse_data$point_of_collaspe[i]<-tempp
# }
# 
# 
# 
# 
# #plotting point of collapse
# tempdat<-point_of_collapse_data #%>% filter(no_of_patches !=20)
# 
# tempdat$point_of_collaspe<-tempdat$point_of_collaspe
# hist(tempdat$point_of_collaspe)
# qqnorm((tempdat$point_of_collaspe)^4)
# 
# summary(model1<-glm(point_of_collaspe/5 ~ no.of.species*no_of_patches+a,
#                     family = quasibinomial,
#                     data=tempdat))
# 
# #rule of thumb of over dispersion: the ration of residual deviance to df should be 1. 
# 
# 
# point_of_collapse_data %>%  
#   filter_all(all_vars(!is.infinite(.))) %>% 
#  # filter(no_of_patches !=20) %>% 
#   filter(type.of.collapse ==  "global") %>% 
#   ggplot(aes(x = no.of.species,
#              y=point_of_collaspe/5, 
#              color = factor(no_of_patches)))+
#   #geom_point(outlier.size = 0)+
#   geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
#   ylab("Point of transition")+
#   xlab("Network size")+
#   stat_smooth(method = "glm",se = T,
#               method.args = list(family = "quasibinomial"))+
#   scale_color_brewer(palette="Dark2")+
#   theme_bw()+
#   facet_grid(.~a)
# 
# 
# 
# #area under the curve for signals of collapse
# 
# auc_ews_data<-expand.grid(`a`= c(0,0.05,0.15),
#                           `no.of.species`=as.numeric(levels((as.factor(fact$network.size)))),
#                           `no_of_patches` =c(2,5,10,20),
#                           `type.of.collapse` =  "global",
#                           `random_seed`=4327+(1:1)*100) %>% 
#   as_tibble %>% 
#   mutate(`auc_regional_variability`=0,
#          `auc_spatial_correlation`=0,
#          `auc_alpha_variation`=0,
#          `auc_spatial_variation`=0,
#          `auc_spatial_sd`=0,
#          `auc_sd`=0,
#          `auc_ar1`=0,
#          `spatial_ar1`=0)
# 
# 
# for(i in 1:nrow(auc_ews_data)){
# dtt<-aucsignal(dat = fact, 
#                mut_strength = fact$g0,
#           type_of_collapse =auc_ews_data$type.of.collapse[i],
#           patches = auc_ews_data$no_of_patches[i],
#           species = auc_ews_data$no.of.species[i],
#             dispersal_rate = auc_ews_data$a[i],
#           seed = auc_ews_data$random_seed[i],
#           tipping.point = point_of_collapse_data$point_of_collaspe[i])
# 
#   
#   auc_ews_data$auc_regional_variability[i] <-dtt$auc_regional_variability
#   auc_ews_data$auc_spatial_correlation[i] <- dtt$auc_spatial_correlation
#   auc_ews_data$auc_spatial_variation[i]<-dtt$auc_spatial_variation
#   auc_ews_data$auc_spatial_sd[i] <-dtt$auc_spatial_sd
#   auc_ews_data$auc_sd[i] <-dtt$auc_sd
#   auc_ews_data$auc_ar1[i] <- dtt$auc_ar1
#   auc_ews_data$auc_alpha_variation[i] <-dtt$auc_alpha_variation
#   auc_ews_data$spatial_ar1[i] <- dtt$auc_spatial_ar1
#   auc_ews_data$web
#    
# }
# 
# # auc ews 
# tempauc<-auc_ews_data #%>% filter(no_of_patches !=20)
# 
# summary(model1<-lm(log(auc_spatial_variation) ~ no.of.species*no_of_patches*a,
#                     data=tempauc))
# hist(log(tempauc$auc_spatial_variation))
# plot(model1)
# 
# auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
#   #filter(no_of_patches != 20) %>% 
#   filter(type.of.collapse ==  "global") %>% 
#   ggplot(aes(x = no.of.species,
#              y=auc_spatial_variation, 
#              color = factor(no_of_patches)))+
#   #geom_boxplot(outlier.size = 0)+
#   geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
#   ylab("AUC spatial variation")+
#   xlab("Network size")+
#   stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
#   scale_color_brewer(palette="Dark2")+
#   theme_bw()+
#   facet_grid(.~a)
# 
# 
# 
# # auc standard deviation
# auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
#   #filter(no_of_patches != 20) %>% 
#   filter(type.of.collapse ==  "global") %>% 
#   ggplot(aes(x = no.of.species,
#              y=auc_sd, 
#              color = factor(no_of_patches)))+
#   #geom_boxplot(outlier.size = 0)+
#   geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
#   ylab("AUC species level standard deviation ")+
#   xlab("Network size")+
#  # stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
#   scale_color_brewer(palette="Dark2")+
#   theme_bw()+
#   facet_grid(.~a)
# 
# 
# #auc alpha variability
# auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
#   #filter(no_of_patches != 20) %>% 
#   filter(type.of.collapse ==  "global") %>% 
#   ggplot(aes(x = no.of.species,
#              y=auc_alpha_variation, 
#              color = factor(no_of_patches)))+
#   geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
#   ylab("AUC alpha variation")+
#   xlab("Network size")+
#   #stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
#   scale_color_brewer(palette="Dark2")+
#   theme_bw()+
#   facet_grid(.~a)
# #auc_regional_variability
# 
# 
# auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
#   filter(type.of.collapse ==  "global") %>% 
#   #filter(no_of_patches != 20) %>% 
#   ggplot(aes(x = no.of.species,
#              y=auc_regional_variability, 
#              color = factor(no_of_patches)))+
#   geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
#   ylab("AUC regional variation")+
#   xlab("Network size")+
#   stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
#   scale_color_brewer(palette="Dark2")+
#   theme_bw()+
#   facet_grid(.~a)
# 
# #spatial autocorrelation at first lag1
# auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
#   ggplot(aes(x = no.of.species,
#              y=spatial_ar1, 
#              color = factor(no_of_patches)))+
#   geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
#   ylab("AUC regional variation")+
#   xlab("Network size")+
#   #stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
#   scale_color_brewer(palette="Dark2")+
#   theme_bw()+
#   facet_grid(.~a)
# 
# 
# #spatial autocorrelation at first lag1
# 
# summary(model1<-lm((auc_ar1) ~ no.of.species+no_of_patches+a,
#                    data=tempauc))
# hist(log(tempauc$auc_spatial_variation))
# plot(model1)
# 
# 
# auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
#   ggplot(aes(x = no.of.species,
#              y=auc_ar1, 
#              color = factor(no_of_patches)))+
#   geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
#   ylab("AUC regional variation")+
#   xlab("Network size")+
#   stat_smooth(method = "lm", se=F, formula = y ~x)+
#   scale_color_brewer(palette="Dark2")+
#   theme_bw()+
#   facet_grid(.~a)
# 
# 
# #spatial autocorrelation at first lag1
# auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
#   ggplot(aes(x = no.of.species,
#              y=auc_spatial_correlation, 
#              color = factor(no_of_patches)))+
#   geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
#   ylab("AUC regional variation")+
#   xlab("Network size")+
#   #stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
#   scale_color_brewer(palette="Dark2")+
#   theme_bw()+
#   facet_grid(.~a)
# 
# #creating network size classes 
# 
