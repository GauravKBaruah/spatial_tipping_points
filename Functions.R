
adj.mat<-function(data){
  #dat <- paste('network.csv',sep='')
  d <- read.csv(file=data,header=FALSE )
  dat<-as.matrix(d)
  dat[dat > 0] = 1
  dat<-apply(dat,2,as.numeric)
  return(dat)}


grid_arrange_shared_legend <- function(..., ncol, nrow, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol =ncol , nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}

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


Connectance<-function(web)
{
  return(sum(web)/(ncol(web)*nrow(web)))}



collapse.point<-function(dat,mut_strength,type_of_collapse, patches, species ,dispersal_rate,
                         seed){
  
  
  temp<-dat %>% filter(type.of.collapse == type_of_collapse,
                       no_of_patches== patches, 
                       network.size==species, 
                       a==dispersal_rate,
                       random_seed==seed,
                       metacomm_biomass==metacomm_biomass)
  
  
  
  rate.of.change.biomass<-abs(temp$metacomm_biomass[1:length(temp$metacomm_biomass)-1] -  
                                temp$metacomm_biomass[2:length(temp$metacomm_biomass)])
  rate.of.change.strength<- 0.15
  nestedness<-temp$nestedness[1]
  conn<-temp$connectance[1]
  
  
  
  tipping.point.index<- which(temp$metacomm_biomass < 0.2*temp$metacomm_biomass[1])[1]
    

  tipping.point <-  (temp$g0[tipping.point.index])
  
  point.of.collapse<- tipping.point
  if(point.of.collapse == -Inf ){ point.of.collapse <- NA}
  
  
  return(list(point.of.collapse=point.of.collapse,nestedness =nestedness , connectance = conn ))
  
}


aucsignal<-function(dat,mut_strength,type_of_collapse, patches, species ,dispersal_rate,
                    seed){
  
  temp<-dat %>% filter(type.of.collapse == type_of_collapse,
                       no_of_patches== patches, 
                       network.size==species, 
                       a==dispersal_rate,
                       random_seed==seed)
  
  nestedness<-temp$nestedness[1]
  connectance<-temp$connectance[1]
  
  auc_ar1<-auc_sd<-auc_spatial_sd<-auc_spatial_variation<-auc_spatial_ar1<-auc_alpha_variation<-auc_spatial_correlation<-auc_regional_variability<-numeric()
  
  
  #temp<-te %>% filter(random_seed ==unique(te$random_seed)[i] )
  # temp<-na.omit(temp)
  #if(tipping.point == 0){
    index_tp<-length(temp$g0)
  #}else { index_tp <- length(temp$g0)
  #}
  #area under the curve for species level ar1
  auc_ar1<-abs(AUC(temp$g0[1:index_tp],(temp$species.ar1[1:index_tp])))
  
  #area under the curve for spatial autocorrelation ar1
  auc_spatial_ar1<-abs(AUC(temp$g0[1:index_tp],(temp$spatial.ar1[1:index_tp])))
  
  
  #area under the curve for species level sd
  auc_sd<-abs(AUC(temp$g0[1:index_tp], temp$species.sd[1:index_tp]))
  
  #area under the curve for spatial  sd
  auc_spatial_sd<-abs(AUC(temp$g0[1:index_tp], temp$spatial.sd[1:index_tp]))
  
  #area under the curve for spatial variability
  auc_spatial_variation<-abs(AUC(temp$g0[1:index_tp], temp$beta_variability[1:index_tp]))
  
  #area under the curve for alpha variaibility
  auc_alpha_variation<-abs(AUC(temp$g0[1:index_tp], temp$alpha_variability[1:index_tp]))
  
  #area under the curve for local level synchrony
  auc_alpha_variation<-abs(AUC(temp$g0[1:index_tp], temp$local_synchrony[1:index_tp]))
  
  #area under the curve for spatial correlation
  auc_spatial_correlation<-abs(AUC(temp$g0[1:index_tp], temp$spatial.correlation[1:index_tp]))
  
  #area under the curve for regional variability
  auc_regional_variability<-abs(as.numeric(AUC(temp$g0[1:index_tp], temp$metacommunity_variability[1:index_tp])))
  
  
  
  return(list(auc_regional_variability=auc_regional_variability,
              auc_spatial_correlation=auc_spatial_correlation,
              auc_alpha_variation=auc_alpha_variation,
              auc_spatial_variation=auc_spatial_variation,
              auc_spatial_sd=auc_spatial_sd,
              auc_spatial_ar1=auc_spatial_ar1,
              auc_sd=auc_sd,
              auc_ar1=auc_ar1,
              nestedness=nestedness,
              connectance=connectance))
  
}




Morans.i<-function(t,state, M, no.of.patches){
  dtemp<-numeric()
  patch.no<-length(state[1,])
  W= patch.no-1 
  temp<-matrix(NA,patch.no,patch.no)
  for (i in 1:no.of.patches){
    for (j in 1: (no.of.patches)){
      
      temp[i,j]<- M[i,j]*(state[t,i] - mean(state[t,]))*(state[t,j] - mean(state[t,])) # / ( (state[t,i]-mean(state[t,]) )^2)
      
    }
    dtemp[i]<- sum(temp[i,])/(1/patch.no * sum( (state[t,] - mean(state[t,]) )^2 ))
    
  }
  
  
  
  I_t <- (1/W)*sum(dtemp) 
  return(I_t)
}


#main dynamical funciton
#time: total time,
#state: abundances 
#pars: list of parameters


eqs<-function(time, state ,pars){
  
  Aspecies<-pars$Aspecies
  Pspecies<-pars$Pspecies
  no.of.patches<-pars$patches
  Disp.mat<-pars$Disp.mat
  Ga<-pars$Ga
  Gp<-pars$Gp
  alpha_p<-pars$alpha_p
  alpha_a<-pars$alpha_a
  gamma<-pars$gamma
  ra<-pars$ra
  rp<-pars$rp
  dt<-0.05
  N<-array(NA,dim=c(time,no.of.patches,Aspecies))
  Np<-array(NA,dim=c(time,no.of.patches,Pspecies))
  N[1,,]<-state
  Np[1,,]<-state
  if(pars$type_of_collapse == "global"){
  

  for (t in 1:(time-1)){
    
    for(k in 1:no.of.patches){
      for(i in 1:Aspecies){
        N[t+1,k,i] <- N[t,k,i]+N[t,k,i]*(ra[k,i]- sum(alpha_a[k,i,]*N[t,k,]) + sum(gamma[k,i,]*Np[t,k,]/(1+0.2*gamma[k,i,]*Np[t,k,]) ) )*dt + 
          sum(N[t,,i]*Disp.mat[,k]*pars$Ga - N[t,k,i]*pars$Ga)*dt +rnorm(1,0,0.1)*N[t,k,i]*dt
       }
      
      for(j in 1:Pspecies){
          
        Np[t+1,k,j] <- Np[t,k,j]+Np[t,k,j]*(rp[k,j]- sum(alpha_p[k,j,]*Np[t,k,]) + sum(gamma[k,,j]*N[t,k,]/(1+0.2*gamma[k,,j]*N[t,k,]) ) )*dt + 
          sum(Np[t,,j]*Disp.mat[,k]*pars$Gp - Np[t,k,j]*pars$Gp)*dt + rnorm(1,0,0.1)*dt*Np[t,k,j]
      }
      N[t+1,k,which(N[t+1,k,] < 1e-7)] <- 0
      Np[t+1,k,which(Np[t+1,k,] < 1e-7)] <- 0
      
      }
    }
  }else if(pars$type_of_collapse == "local"){
    
    for (t in 1:(time-1)){
      
      for(k in 1:no.of.patches){
          for(i in 1:Aspecies){
            N[t+1,k,i] <- N[t,k,i]+N[t,k,i]*(ra[k,i]- sum(alpha_a[k,i,]*N[t,k,]) + sum(gamma[k,i,]*Np[t,k,]/(1+0.2*gamma[k,i,]*Np[t,k,]) ) )*dt + 
              sum(N[t,,i]*Disp.mat[,k]*pars$Ga - N[t,k,i]*pars$Ga)*dt +rnorm(1,0,0.25)*dt*N[t,k,i]
            
          }
          
          for(j in 1:Pspecies){
            Np[t+1,k,j] <- Np[t,k,j]+Np[t,k,j]*(rp[k,j]- sum(alpha_p[k,j,]*Np[t,k,]) + sum(gamma[k,j,]*N[t,k,]/(1+0.2*gamma[k,j,]*N[t,k,]) ) )*dt + 
              sum(Np[t,,j]*Disp.mat[,k]*pars$Gp - Np[t,k,j]*pars$Gp)*dt +rnorm(1,0,0.25)*Np[t,k,j]*dt
          }
          N[t+1,k,which(N[t+1,k,]< 1e-6)] <- 0
          Np[t+1,k,which(Np[t+1,k,]< 1e-6)] <- 0
          
          
        } 
        
      
        
    }
    
  }
  
  output= list(N=N,Np=Np)
  return(output)
}


Mcommunity_1 = function(iter, time, ...){
  set.seed(rnorm(1,as.numeric(Sys.time())-Sys.getpid(),10000)) 
  init = time
  replicate =try(eqs(time=init, ...))
  replicate$start = init
  replicate$iter = iter
  return(replicate)
}



#dynamical function that evaluates tipping point at the level of species
eqs_species<-function(time, state ,pars){
  
  Aspecies<-pars$Aspecies
  Pspecies<-pars$Pspecies
  no.of.patches<-pars$patches
  Disp.mat<-pars$Disp.mat
  Ga<-pars$Ga
  Gp<-pars$Gp
  alpha_p<-pars$alpha_p
  alpha_a<-pars$alpha_a
  gamma<-pars$gamma
  ra<-pars$ra
  rp<-pars$rp
  dt<-0.05
  degree.animals<-degree.plants<-numeric()
  N<-array(NA,dim=c(time,no.of.patches,Aspecies))
  Np<-array(NA,dim=c(time,no.of.patches,Pspecies))
  Na_tipping.point<-array(NA,dim=c(no.of.patches,Aspecies))
  Np_tipping.point<-array(NA,dim=c(no.of.patches,Pspecies))
  N_acf_P<-N_avg_P<-array(NA,dim=c(length(pars$g0),no.of.patches,Pspecies))
  N_acf_A<-N_avg_A<-array(NA,dim=c(length(pars$g0),no.of.patches,Aspecies))
  N_sd_P<-array(NA,dim=c(length(pars$g0),no.of.patches,Pspecies))
  N_sd_A<-array(NA,dim=c(length(pars$g0),no.of.patches,Aspecies))
  
  auc_ar1_A<-auc_sd_A<-array(NA, dim=c(no.of.patches,Aspecies))
  auc_ar1_P<-auc_sd_P<-array(NA, dim=c(no.of.patches,Pspecies))
  N[1,,]<-state
  Np[1,,]<-state
  #pars$g0<-c(0.1,2,3,5)
  for(r in 1:length(pars$g0)){
    for (t in 1:(time-1)){
      
      for(k in 1:no.of.patches){
      for(i in 1:Aspecies){
        for(m in 1:Pspecies){
          gamma[k,i,m] <- rnorm(1, pars$g0[r], 0.000)*pars$g[i,m]
        }
      }
      }
      for(k in 1:no.of.patches){
        for(i in 1:Aspecies){
          N[t+1,k,i] <- N[t,k,i]+N[t,k,i]*(ra[k,i]- sum(alpha_a[k,i,]*N[t,k,]) + sum(gamma[k,i,]*Np[t,k,]/(1+0.2*gamma[k,i,]*Np[t,k,]) ) )*dt + 
            sum(N[t,,i]*Disp.mat[,k]*pars$Ga - N[t,k,i]*pars$Ga)*dt +rnorm(1,0,0.05)*N[t,k,i]*dt
        }
        
        for(j in 1:Pspecies){
          
          Np[t+1,k,j] <- Np[t,k,j]+Np[t,k,j]*(rp[k,j]- sum(alpha_p[k,j,]*Np[t,k,]) + sum(gamma[k,,j]*N[t,k,]/(1+0.2*gamma[k,,j]*N[t,k,]) ) )*dt + 
            sum(Np[t,,j]*Disp.mat[,k]*pars$Gp - Np[t,k,j]*pars$Gp)*dt + rnorm(1,0,0.05)*dt*Np[t,k,j]
        }
        N[t+1,k,which(N[t+1,k,] < 1e-6)] <- 0
        Np[t+1,k,which(Np[t+1,k,] < 1e-6)] <- 0
        
      }
    }
    
    for(k in 1:no.of.patches){
      for(j in 1:Pspecies){
        N_acf_P[r,k,j]<-acf(Np[100:time,k,j], lag.max = 1, type = c("correlation"), plot = FALSE)$acf[2]
        N_sd_P[r,k,j]<-sd(Np[100:time,k,j],na.rm=T)
        N_avg_P[r,k,j]<-mean(Np[100:time,k,j])
      }
      
      for(l in 1:Aspecies){
        N_acf_A[r,k,l]<- acf(N[100:time,k,l], lag.max = 1, type = c("correlation"), plot = FALSE)$acf[2]
        N_sd_A[r,k,l]<- sd(N[100:time,k,l],na.rm=T)
        
          N_avg_A[r,k,l]<-mean(N[100:time,k,l])
      }
  }
  

  }
  
  
  
  for(k in 1:no.of.patches){
    for(j in 1:Pspecies){
    rate_of_change<- abs(N_avg_P[1:length(pars$g0),k,j] -  N_avg_P[2:(length(pars$g0)-1),k,j])
    tipping.point.index<- pars$g0[which(N_avg_P[1:length(pars$g0),k,j] < 0.2*N_avg_P[1,k,j])][1]
    tipping.point.index[!is.finite(tipping.point.index)] <- NA
    Np_tipping.point[k,j] <- tipping.point.index
   # P_Acf[k,j]<- acf(Np[100:time,k,j], lag.max = 1, type = c("correlation"), plot = FALSE)$acf[2]
  #  P_sd[k,j] <- sd(Np[100:time,k,j],na.rm=T)
   

    auc_ar1_P[k,j]<-abs(AUC(pars$g0,(N_acf_P[,k,j])))
    auc_sd_P[k,j]<-abs(AUC(pars$g0,(N_sd_P[,k,j])))
    degree.plants[j]<-sum(pars$g[,j])
    }
    for(l in 1:Aspecies){
      rate_of_change<- abs(N_avg_A[1:length(pars$g0),k,l] -  N_avg_A[2:(length(pars$g0)-1),k,l])
      tipping.point.index<-  pars$g0[which(N_avg_A[1:length(pars$g0),k,l] < 0.2*N_avg_A[1,k,l])][1]
      tipping.point.index[!is.finite(tipping.point.index)] <- NA
      Na_tipping.point[k,l] <- tipping.point.index
     # A_Acf[k,l]<- acf(N[100:time,k,l], lag.max = 1, type = c("correlation"), plot = FALSE)$acf[2]
    #  A_sd[k,l] <- sd(N[100:time,k,l],na.rm=T)
      auc_ar1_A[k,l]<-abs(AUC(pars$g0,(N_acf_A[,k,l])))
      auc_sd_A[k,l]<-abs(AUC(pars$g0,(N_sd_A[,k,l])))
      degree.animals[l]<-sum(pars$g[l,])
    }
    
  }
  Np_tipping_point_1=colMeans(Np_tipping.point,na.rm=T)
  Na_tipping.point_2=colMeans(Na_tipping.point,na.rm=T)
  A_acf<-colMeans(auc_ar1_A,na.rm=T)
  P_acf<-colMeans(auc_ar1_P,na.rm = T)
  A_sd<-colMeans(auc_sd_A,na.rm=T)
  P_sd<-colMeans(auc_sd_P,na.rm=T)
  
  return(list(Np_tipping_point_1=Np_tipping_point_1,
              Na_tipping.point_2=Na_tipping.point_2, 
              degree.animals=degree.animals,
              degree.plants=degree.plants,
              A_acf=A_acf,
              P_acf=P_acf,
              A_sd=A_sd,
              P_sd=P_sd))
  
  
}


Mcommunity_sp = function(iter, time, ...){
  set.seed(rnorm(1,as.numeric(Sys.time())-Sys.getpid(),10000)) 
  init = time
  replicate =try(eqs_species(time=init, ...))
  replicate$start = init
  replicate$iter = iter
  return(replicate)
}



#function that estiamtes all the statistical metrics at the community level, metacommunity level and species level such as 
# alpha variability, beta variability, average species autocorrelation, standard deviation.

meta.stats<-function(dat, patches, species, time,disp.mat, type_of_collapse,patch_affected){
  
if(type_of_collapse == "global"){
   data<-array(NA,dim=c(time,patches,species))
  for(k in 1:patches){
      data[,k,]<-cbind(dat[[1]]$N[,k,], dat[[1]]$Np[,k,])
  }
 }
  a<-b<-g<-numeric()
  temp<-matrix(NA,patches,patches)
  v1<-array(NA,dim=c(species,species,patches,patches))
  v<-array(NA,dim=c(species,species,patches,patches))
  t1<-d1<-array(NA,dim=c(species,species,patches))
  t2<-d2<-array(NA,dim=c(species,species))
  biomass.patches<-array(NA,dim=c(time,patches))
  t3<-d3<-numeric()
  for(t in 1:time){
    for (i in 1:patches){
      a[i] <- length(which(data[t,i,]>0))
      biomass.patches[t,i] <- sum(data[t,i,])
    }
    b[t]<-sum(a,na.rm=T)/patches
    g[t]<-sum(a,na.rm=T)
    }
  temp<-mat<-species_Acf<-species_sd<-array(NA,dim=c(patches,species))
  u<-cv_patches<-biomass_in_patches<-ACF_patches<-biomass_patches<-numeric()
  for( i in 1:patches){
    for(j in 1:species){
      temp[i,j] <- sum(data[,i,j])/time
      mat[i,j] <- data[time,i,j]
      species_Acf[i,j]<- acf(data[100:time,i,j], lag.max = 1, type = c("correlation"), plot = FALSE)$acf[2]
      species_sd[i,j] <- sd(data[100:time,i,j],na.rm=T)
      
    }
    u[i]<-sum(temp[i,])
    cv_patches[i]<-  sd(biomass.patches[100:time,i])/mean(biomass.patches[100:time,i])
    ACF_patches[i]<- acf(biomass.patches[100:time,i], lag.max = 1, type = c("correlation"), plot = FALSE)$acf[2]
    biomass_in_patches[i] <- sum(temp[i,],na.rm=T)
    biomass_patches[i]<- mean(biomass.patches[1:time,i],na.rm=T)
  }
  
  #alpha variability
  avg.cv_patch <-sum(cv_patches,na.rm=T)/sum(biomass_in_patches,na.rm=T)
  
  #community autocorrelation
  avg.acf_patch<-mean(ACF_patches,na.rm=T)
  
  temp2<-matrix(NA, patches, patches)
  species_level_acf<-species_level_sd<-dtemp<-disp_row<-numeric()
  diag(disp.mat)<-0
  for(i in 1:species){
    for(j in 1:species){
      for(k in 1:patches){
        for(l in 1:patches){
          v[i,j,k,l] <- sum( (data[,k,i]- temp[k,i])*(data[,l,j] - temp[l,j]))/(time-1)
          v1[i,i,k,k] <- sum( sqrt((data[,k,i]- temp[k,i])*(data[,k,i] - temp[k,i])))/(time-1)
          temp2[k,l]<- disp.mat[k,l]*(biomass.patches[time,k] - mean(biomass.patches[time,]))*(biomass.patches[time,l] - mean(biomass.patches[time,])) # / ( (state[t,i]-mean(state[t,]) )^2)
        }
        
        dtemp[k]<- patches*sum(temp2[k,])/(sum( (biomass.patches[time,] - mean(biomass.patches[time,]) )^2 ))
        disp_row[k] <- sum(disp.mat[k,])
        
        
        t1[i,j,k]<-sum(v[i,j,k,],na.rm=T)
        d1[i,i,k]<-sum(v[i,i,k,],na.rm = T)
      }
      t2[i,j]<-sum(t1[i,j,])
     
    }
    #d3[i]<-sum(d2[i,])
    t3[i]<-sum(t2[i,])
    d2[i,i]<-sum(d1[i,i,],na.rm = T)
    species_level_acf[i]<-mean(species_Acf[,i],na.rm=T)
    species_level_sd[i]<-mean(species_sd[,i],na.rm=T)
  }
  
  I_t <-  sum(dtemp,na.rm=T)/sum(disp_row,na.rm=T)   # (1/(patches-1))*sum(dtemp) 
  
  viikk<-vijkk<-array(NA,dim=c(patches,patches))
  v_sig_kk<-array(NA,dim=c(species,patches,patches))
  for(k in 1:patches){
    
      viikk[k,k]<- sum(sqrt(v1[,,k,k]), na.rm=T)
    for(i in 1:species){
      v_sig_kk[i,k,k]<- sum(v[i,,k,k],na.rm = T)
    }
      vijkk[k,k]<-sum(v_sig_kk[,k,k],na.rm = T)
  }
  
  avg.local.synchrony<-sum(vijkk,na.rm = T)/sum(viikk,na.rm = T)
  metacom_variability<-sqrt(sum(t3))/sum(u)
  reg.synchrony<-metacom_variability/sum(d2,na.rm = T)
  temporal.gamma.richness<-b[time]
  alpha.richness<-sum(b,na.rm=T)/time
  gamma.richness<-sum(g,na.rm=T)/time
  spatial.correlation <-I_t
  spatial.beta.richness<- gamma.richness-alpha.richness
  alpha.variability <- avg.cv_patch
  alpha_acf<-avg.acf_patch
  if(type_of_collapse== "global"){
  mean.biomass.patch.affected<-sum(biomass_patches)
  }else if(type_of_collapse == "local"){
    mean.biomass.patch.affected<-mean.patch.affected
  }
  beta.variability <- avg.cv_patch/metacom_variability
  #compositional_dissimilarity<-vegdist(mat,method="bray")
  metacomm_biomass<-sum(biomass_patches)
    
  return(list(temporal.gamma.richness=temporal.gamma.richness,
         alpha.richness=alpha.richness,
         mean.biomass.patch.affected=mean.biomass.patch.affected,
         gamma.richness=gamma.richness,
         spatial.correlation=spatial.correlation,
         Acf_patch=alpha_acf, #temporal autocorrelation across patches
         species_level_acf=species_level_acf, #autocorrelation at species level
         species_level_sd=species_level_sd, #cv at species level
         spatial.beta.richness=spatial.beta.richness,
         metacom_variability=metacom_variability, #temporal variability at the metacommunity level
         alpha.variability=alpha.variability, #temporal variability at the patch level
         beta.variability=beta.variability, #spatial variability
         #compositional_dissimilarity=compositional_dissimilarity,
         metacomm_biomass=metacomm_biomass,
         avg.local.synchrony=avg.local.synchrony,
         reg.synchrony=reg.synchrony
         ))
  }

