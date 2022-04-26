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


grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  # grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}
#1. Web of interaction as matrix
#2. species trait variances

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


mydir = 'plant_pollinator'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
myfiles<-myfiles

g0 = seq(5,0,-0.2)

#creating an empty data frame


fact<- expand.grid(`a`= c(0,0.05,0.15),
                   `web`=myfiles[1:53],
                   `no_of_patches` =c(2,5,10,20),
                   `type.of.collapse` =  "global",
                   `random_seed`=4327+(1:1)*100) %>%
  as_tibble %>%
  mutate(`tipping.point`=0,
         degree.plants=0,
         degree.animals=0,
         nestedness=0,
         connectance=0,
         network.size=0)

for(i in 1:nrow(fact)){
  
  g<-adj.mat(myfiles[which(myfiles == fact$web[i])]) #network web names
  #g<-g[-1,-1] 
  fact$network.size[i]<- nrow(g)+ncol(g)
  fact$nestedness[i] <- nestedness_NODF(g)
  fact$connectance[i]<-Connectance(g)
  
}

fact<-fact %>% filter(web != "plant_pollinator/M_PL_046.csv" & web != "plant_pollinator/M_PL_013.csv" )


new_df<-NULL

for(j in 423:nrow(fact)){
  
  if(fact$no_of_patches[j] == 2){
    coordinates<- Cordinates_list$coordinates_2
  }else if(fact$no_of_patches[j] == 5){
    coordinates<- Cordinates_list$coordinates_5
  }else if(fact$no_of_patches[j] == 10){
    coordinates<- Cordinates_list$coordinates_10
  }else if(fact$no_of_patches[j] == 20){
    coordinates<- Cordinates_list$coordinates_20
  }
  
  
  # if(fact$no_of_patches[j] == 2 & fact$type.of.collapse[j] == "local"){
  #   coordinates<- Cordinates_list$coordinates_2
  #   patch_affected_2<-patch_removed_2
  # }else if(fact$no_of_patches[j] == 5 & fact$type.of.collapse[j] == "local"){
  #   coordinates<- Cordinates_list$coordinates_5
  #   patch_affected_5<- patch_removed_5
  # }else if(fact$no_of_patches[j] == 10 & fact$type.of.collapse[j] == "local"){
  #   coordinates<- Cordinates_list$coordinates_10
  #   patch_affected_15 <- patch_removed_15
  # }else if(fact$no_of_patches[j] == 20 & fact$type.of.collapse[j] == "local"){
  #   coordinates<- Cordinates_list$coordinates_20
  #   patch_affected_20<- patch_removed_20
  # }
  
  
  Disp.mat<-BuildRndLandscape(N = fact$no_of_patches[j],w=0.25,Coordinates=coordinates,
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
  gamma <-array(NA,dim=c(fact$no_of_patches[j], Aspecies, Pspecies ))
  
  
  
  for(k in 1:fact$no_of_patches[j]){
    alpha_a[k,,] <- runif(Aspecies^2, 0.005,0.01)
    alpha_p[k,,] <-runif(Pspecies^2,0.005,0.01)
    diag(alpha_a[k,,])<-0.5;diag(alpha_p[k,,])<-0.5
    
    # }else{ alpha_a <-alpha_a;alpha_p<-alpha_p
    # }
    # 
    
    ra[k,] <- runif(Aspecies,-0.1,-0.01)
    rp[k,] <- runif(Pspecies,-0.1,-0.01)

  
  
  Ga <- fact$a[j]
  Gp <- fact$a[j]
  
  }
  pars<-list(Aspecies=Aspecies,Pspecies=Pspecies,patches=fact$no_of_patches[j], Disp.mat=Disp.mat,Ga=Ga,Gp=Gp,g0=g0,
             alpha_a=alpha_a,alpha_p=alpha_p,ra=ra,rp=rp,g=g,gamma=gamma,model="Thompson_et_al",type_of_collapse= fact$type.of.collapse[j])
  
  start.time = 700
  state<-0.5
  model.t<-lapply(1, Mcommunity_sp,time=start.time,state=state,
                  pars=pars)
  
 
  webname<-fact$web[j]
  
  ddf <-as.data.frame(cbind( rep(as.character(webname),each=(nrow(g)+ncol(g))), 
                       rep(fact$a[j],each=(nrow(g)+ncol(g))),
                       rep(fact$no_of_patches[j],each=(nrow(g)+ncol(g))),
                       rep(seq(1,(nrow(g)+ncol(g)),1)),
                       c(model.t[[1]]$Na_tipping.point_2,model.t[[1]]$Np_tipping_point_1),
                       c(model.t[[1]]$degree.animals,model.t[[1]]$degree.plants),
                       c(model.t[[1]]$A_acf,model.t[[1]]$P_acf),
                       c(model.t[[1]]$A_sd,model.t[[1]]$P_sd),
                       rep(nestedness_NODF(g), each=((nrow(g)+ncol(g))) ),
                       rep(Connectance(g), each=((nrow(g)+ncol(g))) )))
    colnames(ddf)<-c("web","a","patches","species","tipping.points","degree","Autocorrelation","sd","nestedness", "connectance")
    
    new_df<-rbind(new_df, ddf )
   
  
  
  print(j)
  }
  
 save(new_df,file="Species_tipping_points_spatial_data_apr15.RData")


#15019
#
#save(fact, file="scale_of_ews_real_webs_20may.RData")
#load("scale_of_ews_real_webs_20may.RData")
#load("scale_of_ews_real_webs_15may.RData")
#fact<-rbind(fact,fact_1)

#save(fact,file="scale_of_ews_webs_all_data.RData")

load("scale_of_ews_webs_all_data.RData")
fact<-fact %>% filter(web != "plant_pollinator/M_PL_013.csv" & web!= "plant_pollinator/M_PL_046.csv")

#point of collapse estimation
point_of_collapse_data<-expand.grid(`a`= c(0,0.05,0.15),
                                    `no.of.species`=as.numeric(levels(as.factor(fact$network.size))),
                                    `no_of_patches` =c(2,5,10,20),
                                    `type.of.collapse` = "global",
                                    `random_seed`=4327+(1:1)*100) %>% 
  as_tibble %>% 
  mutate(`point_of_collaspe`=0,
         nestedness=0,
         connectance=0)

for(i in 1:nrow(point_of_collapse_data)){
  tempp<-collapse.point(dat = fact, 
                        mut_strength = fact$g0, 
                        type_of_collapse = point_of_collapse_data$type.of.collapse[i], 
                        patches = point_of_collapse_data$no_of_patches[i],
                        species = point_of_collapse_data$no.of.species[i], 
                        dispersal_rate = point_of_collapse_data$a[i] , 
                        seed = point_of_collapse_data$random_seed[i])
  
  point_of_collapse_data$point_of_collaspe[i]<-tempp$point.of.collapse
  point_of_collapse_data$nestedness[i]<-tempp$nestedness
  point_of_collapse_data$connectance[i]<-tempp$connectance
}




#plotting point of collapse
tempdat<-point_of_collapse_data #%>% filter(no_of_patches !=20)

tempdat$point_of_collaspe<-tempdat$point_of_collaspe
hist(tempdat$point_of_collaspe)
qqnorm((tempdat$point_of_collaspe)^4)

summary(model1<-glm(point_of_collaspe/5 ~ no.of.species*no_of_patches+a,
                    family = quasibinomial,
                    data=tempdat))

#rule of thumb of over dispersion: the ration of residual deviance to df should be 1. 


point_of_collapse_data %>%  
  filter_all(all_vars(!is.infinite(.))) %>% 
  # filter(no_of_patches !=20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = nestedness,
             y=log(point_of_collaspe+1), 
             color = factor(no_of_patches)))+
  #geom_point(outlier.size = 0)+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  xlab("Network size")+
  ylab( expression(paste("log(threshold mutualistic strength),"," ", gamma[0])))+
 stat_smooth(method = "lm",se = T,
              formula = y ~ x)+
  scale_color_brewer(palette="Dark2")+
  labs(color="No. of patches")+
  theme_bw()+
  facet_grid(.~a)


point_of_collapse_data %>%  
  filter_all(all_vars(!is.infinite(.))) %>% 
  # filter(no_of_patches !=20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = connectance,
             y=log(point_of_collaspe+1), 
             color = factor(no_of_patches)))+
  #geom_point(outlier.size = 0)+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  xlab("Nestedness")+
  ylab( expression(paste("threshold mutualistic strength,"," ", gamma[0])))+
  #stat_smooth(method = "lm",se = F,
   #           formula = y ~ x )+
  scale_color_brewer(palette="Dark2")+
  labs(color="No. of patches")+
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
         nestedness=0,
         connectance=0)


for(i in 1:nrow(auc_ews_data)){
  dtt<-aucsignal(dat = fact, 
                 mut_strength = fact$g0,
                 type_of_collapse =auc_ews_data$type.of.collapse[i],
                 patches = auc_ews_data$no_of_patches[i],
                 species = auc_ews_data$no.of.species[i],
                 dispersal_rate = auc_ews_data$a[i],
                 seed = auc_ews_data$random_seed[i])
  
  
  auc_ews_data$auc_regional_variability[i] <-dtt$auc_regional_variability
  auc_ews_data$auc_spatial_correlation[i] <- dtt$auc_spatial_correlation
  auc_ews_data$auc_spatial_variation[i]<-dtt$auc_spatial_variation
  auc_ews_data$auc_spatial_sd[i] <-dtt$auc_spatial_sd
  auc_ews_data$auc_sd[i] <-dtt$auc_sd
  auc_ews_data$auc_ar1[i] <- dtt$auc_ar1
  auc_ews_data$auc_alpha_variation[i] <-dtt$auc_alpha_variation
  auc_ews_data$spatial_ar1[i] <- dtt$auc_spatial_ar1
  auc_ews_data$connectance[i]<-dtt$connectance
  auc_ews_data$nestedness[i]<-dtt$nestedness
  
}

# auc ews 
tempauc<-auc_ews_data #%>% filter(no_of_patches !=20)

summary(model1<-lm(log(auc_spatial_variation) ~ connectance+ no_of_patches+a,
                   data=tempauc))
hist(log(tempauc$auc_spatial_variation))
plot(model1)

auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  #filter(no_of_patches != 20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = no.of.species,
             y=log(auc_spatial_variation), 
             color = factor(no_of_patches)))+
  #geom_boxplot(outlier.size = 0)+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("log(AUC spatial variation)")+
  xlab("Network size")+
  stat_smooth(method = "lm", se=T)+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  labs(color="No. of patches")+
  facet_grid(.~a)



auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  #filter(no_of_patches != 20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = connectance,
             y=log(auc_spatial_variation), 
             color = factor(no_of_patches)))+
  #geom_boxplot(outlier.size = 0)+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("log(AUC spatial variation)")+
  xlab("Network connectance")+
  stat_smooth(method = "lm", se=T)+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  labs(color="No. of patches")+
  facet_grid(.~a)



# auc standard deviation


summary(model1<-lm((auc_sd) ~ nestedness+ no_of_patches+a,
                   data=tempauc))
#hist(log(tempauc$auc_spatial_variation))
plot(model1)

auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  #filter(no_of_patches != 20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = connectance,
             y=auc_sd, 
             color = factor(no_of_patches)))+
  #geom_boxplot(outlier.size = 0)+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC species level standard deviation ")+
  xlab("Network size")+
   stat_smooth(method = "lm", se=F,formula = y ~ x)+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  labs(color="No. of patches")+
  facet_grid(.~a)




#auc alpha variability
summary(model1<-lm(log(auc_alpha_variation+1) ~ no.of.species+ factor(no_of_patches)+a,
                   data=auc_ews_data))
hist((tempauc$auc_alpha_variation))
plot(model1)


auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  #filter(no_of_patches != 20) %>% 
  filter(type.of.collapse ==  "global") %>% 
  ggplot(aes(x = no.of.species,
             y=log(auc_alpha_variation+1), 
             color = factor(no_of_patches) ))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("log(AUC alpha variation)")+
  xlab("Network size")+
  stat_smooth(method = "lm", formula = y ~ x)+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  labs(color="No. of patches")+
  facet_grid(.~a)



#auc_regional_variability

summary(model1<-lm((auc_regional_variability) ~ no.of.species+ no_of_patches+a,
                   data=tempauc))
hist((tempauc$auc_alpha_variation))
plot(model1)


auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  filter(type.of.collapse ==  "global") %>% 
  #filter(no_of_patchesno20) %>% 
  ggplot(aes(x = connectance,
             y=log(auc_regional_variability+1), 
             color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("log(AUC regional variation)")+
  xlab("Network size")+
  stat_smooth(method = "lm", se=T)+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  labs(color="No. of patches")+
  facet_grid(.~a)

#spatial autocorrelation at first lag1
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  ggplot(aes(x = connectance,
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
  ylab("AUC species AR1")+
  xlab("Network size")+
  #stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)


#spatial autocorrelation at first lag1
auc_ews_data %>%  filter_all(all_vars(!is.infinite(.))) %>% 
  ggplot(aes(x = no.of.species,
             y=auc_spatial_correlation, 
             color = factor(no_of_patches)))+
  geom_point(alpha = 5/10, size = 1.5, stroke = 1.5,shape = 21)+
  ylab("AUC spatial correlation")+
  xlab("Network size")+
  #stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  facet_grid(.~a)

#creating network size classes 
fact$network.size_class <- 
  findInterval(fact$network.size, seq(10, max(fact$network.size), 10), left.open = TRUE)




#plotting tipping points for three different networks




#biomass web number 1
b1<-fact %>% filter(web == "plant_pollinator/M_PL_003.csv") %>% 
  ggplot(aes(x=g0, y=metacomm_biomass, color=factor(no_of_patches)))+
  geom_point(alpha =0.7, size =2.5)+xlab("")+ylab("metacommunity biomass")+xlab("mutualistic strength")+
  #geom_smooth(method = 'loess',se=F)+
  theme_bw()+
  theme(legend.position = "none") +
  labs("")+
  facet_grid(.~a)

b2<-fact %>% filter(web == "plant_pollinator/M_PL_060_21.csv") %>% 
  ggplot(aes(x=g0, y=metacomm_biomass, color=factor(no_of_patches)))+
  geom_point(alpha =0.7, size =2.5)+xlab("")+
  ylab("metacommunity biomass")+xlab("mutualistic strength")+
  #geom_smooth(method = 'loess',se=F)+
  theme_bw()+
  theme(legend.position = "none") +
  labs(color=" " )+
  facet_grid(.~a)


b3<-(fact %>% filter(web == "plant_pollinator/M_PL_061_14.csv") %>% 
  ggplot(aes(x=g0, y=metacomm_biomass, color=factor(no_of_patches)))+
  geom_point(alpha =0.7, size =2.5)+xlab("")+
    ylab("metacommunity biomass")+xlab("mutualistic strength")+
  #geom_smooth(method = 'loess',se=F)+
  theme_bw()+
  theme(legend.position = "none") +
  facet_grid(.~a))

grid_arrange_shared_legend(b1,b2,b3,nrow=3,ncol=1)

net1<-adj.mat("plant_pollinator/M_PL_003.csv")
net3<-adj.mat("plant_pollinator/M_PL_061_33.csv")
net2<-adj.mat("plant_pollinator/M_PL_061_18.csv")

par(mfrow=(c(3,1)))
web1<-plotweb(net1,
              method="normal",ybig=0.1, y.width.low = 0.1,
              col.interaction="wheat4",
              bor.col.interaction="white", 
              arrow="no",  col.high="lightblue",
              col.low="tomato",labsize=0.1)
web2<-plotweb(net2,
              method="normal",ybig=0.1, y.width.low = 0.1,
              col.interaction="wheat4",
              bor.col.interaction="white", 
              arrow="no",  col.high="lightblue",
              col.low="tomato",labsize=0.1)
web3<-plotweb(net3,
              method="normal",ybig=0.1, y.width.low = 0.1,
              col.interaction="wheat4",
              bor.col.interaction="white",
              arrow="no",  col.high="lightblue",
              col.low="tomato",labsize=0.1)

grid_arrange_shared_legend(plot1,plot2,
                           plot3,plot4,
                           plot5,plot6,
                           nrow=3,ncol=2)

#species richness

fact %>% 
  ggplot(aes(x=g0, y=patch_richness, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("metacommunity biomass")+xlab("mutualistic strength")+
  geom_smooth(method = 'loess',se=TRUE)+theme_bw()+
  facet_grid(a~network.size_class)


#species level ar1
fact %>%  filter(type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=species.ar1, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("species ar1")+
  geom_smooth(method = 'loess',se=TRUE)+theme_bw()+facet_grid(a~network.size_class)


# temporal patch variability
fact %>%  filter(type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=alpha_variability, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("patch variability")+
  geom_smooth(method = 'loess',se=TRUE)+theme_bw()+facet_grid(a~network.size_class)

# spatial Standard deviation

fact %>%  filter( type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=spatial.sd, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("spatial SD")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~network.size_class)




#spatial correlation


fact %>%  filter(type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=spatial.correlation, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("spatial correlation")+
  geom_smooth(method = 'loess',se=TRUE)+theme_bw()+facet_grid(a~network.size_class)



# regional synchrony



fact %>%  filter(type.of.collapse=="local") %>% 
  ggplot(aes(x=g0, y=regional_synchrony, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("Regional synchrony")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~no.of.species)

fact %>%  filter( type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=regional_synchrony, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("Regional synchrony")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~no.of.species)



# local synchrony

fact %>%  filter(type.of.collapse=="local") %>% 
  ggplot(aes(x=g0, y=local_synchrony, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("local synchrony")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~no.of.species)

fact %>%  filter( type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=local_synchrony, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("local synchrony")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~no.of.species)


# network dissimilarity

fact %>%  filter(type.of.collapse=="local") %>% 
  ggplot(aes(x=g0, y=network_dissimilarity, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("metacommunity variability")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~no.of.species)

fact %>%  filter( type.of.collapse=="global") %>% 
  ggplot(aes(x=g0, y=network_dissimilarity, color=factor(no_of_patches)))+
  geom_point()+xlab("")+ylab("metacommunity variability")+
  geom_smooth(method = 'loess',se=TRUE)+theme_classic()+facet_grid(a~no.of.species)

