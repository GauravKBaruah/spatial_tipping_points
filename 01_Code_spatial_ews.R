rm(list=ls())
library(magic)
library(igraph)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(vegan)
library(NetIndices)
library(betalink)
library(igraph)
library(DescTools)
source('~/Dropbox/EAWAG PostDoc/06_scale_of_EWS/spatial_tipping_points/dispKernels.R',echo=F)
source('~/Dropbox/EAWAG PostDoc/06_scale_of_EWS/spatial_tipping_points/buildlandscape.R', echo=F)
source('~/Dropbox/EAWAG PostDoc/06_scale_of_EWS/spatial_tipping_points/Functions.R', echo=F)




#empty arrays for abundance of species
start.time <-5000


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



pars<-list()


mydir = 'plant_pollinator'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
myfiles<-myfiles[1:56]


#creating an empty data frame
fact<- expand.grid(`w`=0.5,
                   `g0` = seq(6,0,-0.175),
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
        alpha_a[k,,] <- runif(Aspecies^2, 0.005,0.01) #interspecies competition
        alpha_p[k,,] <-runif(Pspecies^2,0.005,0.01) #interspecies competition plants
       diag(alpha_a[k,,])<-0.5;diag(alpha_p[k,,])<-0.5 #intraspecific competition
       
      
        ra[k,] <- runif(Aspecies,-0.1,-0.01) #obligate mutualism 
        rp[k,] <- runif(Pspecies,-0.1,-0.01) #obligate mutualism
      if(fact$type.of.collapse[j] == "global"){
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
  
  start.time = 5000
  state<-0.5
  model.t<-lapply(1, Mcommunity_1,time=start.time,state=state,
                  pars=pars)
  #model.t[[1]]$N[800,1,]
  
  
  if(fact$type.of.collapse[j] == "global") {
    fact$patch_affected[j] <- 0
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
   

  print(j)
}

 #save(fact, file="example_network_collapse.RData")
