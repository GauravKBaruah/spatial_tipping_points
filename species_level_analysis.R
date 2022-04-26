rm(list=ls())
library(magic)
library(igraph)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(vegan)
library(MuMIn)
library(lme4)
library(nlme)
library(cowplot)
library(cowplot)
source('~/Dropbox/EAWAG PostDoc/06_scale_of_EWS/dispKernels.R',echo=F)
source('~/Dropbox/EAWAG PostDoc/06_scale_of_EWS/buildlandscape.R', echo=F)
source('~/Dropbox/EAWAG PostDoc/06_scale_of_EWS/Functions.R', echo=F)
source('~/Dropbox/Zurich PhD Research/1_Chapter_1/rolling_GAMs_methods.R', echo=F)

library(NetIndices)
library(betalink)
library(igraph)
library(DescTools)

adj.mat<-function(data){
  #dat <- paste('network.csv',sep='')
  d <- read.csv(file=data,header=FALSE )
  dat<-as.matrix(d)
  dat[dat > 0] = 1
  dat<-apply(dat,2,as.numeric)
  return(dat)}




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




load("Species_tipping_points_spatial_data_apr15.RData")


head(new_df)
str(new_df)


#converting the dataframe to numeric 
new_df$a<-as.numeric(new_df$a)
new_df$patches<-as.numeric(new_df$patches)
new_df$species <- as.numeric(new_df$species)
new_df$tipping.points<-as.numeric(new_df$tipping.points)
new_df$degree<-as.numeric(new_df$degree)
new_df$nestedness<-as.numeric(new_df$nestedness)
new_df$connectance<-as.numeric(new_df$connectance)
new_df$Autocorrelation<-as.numeric(new_df$Autocorrelation)
new_df$sd<-as.numeric(new_df$sd)

adj.mat<-function(data){
  #dat <- paste('network.csv',sep='')
  d <- read.csv(file=data,header=FALSE )
  dat<-as.matrix(d)
  dat[dat > 0] = 1
  dat<-apply(dat,2,as.numeric)
  return(dat)}
mydir = 'plant_pollinator'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:56]


for( r in 1:53) {
  g<-adj.mat(myfiles[which(myfiles == new_df$web[r])]) #network web names
  #g<-g[-1,-1] 
  new_df$No.species[r]<-nrow(g)+ncol(g)  
  
}
#removing nans
new_df<-na.omit(new_df)



# linear mixed model analysis

#checking for whether linear model is possible
require(car)
require(MASS)
new_df$tipping.points<- new_df$tipping.points+1
qqp(new_df$tipping.points, "lnorm")

hist(log(new_df$tipping.points+1))

#new_df<-new_df %>% filter(web !="plant_pollinator/M_PL_046.csv")
#new_df<-new_df %>% filter(web !="plant_pollinator/M_PL_013.csv")

(model<-lme(log(tipping.points) ~ 
              degree + nestedness  + as.factor(a) , random= ~1|web, data = new_df))

summary(model)

colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
grid.arrange(plot(model,type=c("p","smooth")),
             plot(model,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             qqnorm(model,abline=c(0,1)),ncol=3)

plot(model)
summary(model)

summary(mod<-rlm((tipping.points) ~ 
           degree + nestedness  + a,  method= "MM",data = new_df))


anova(mod)

summary(mod<-lmrob(tipping.points ~ 
                   degree + nestedness  + a,  data = new_df))


#anova(mod)

tdat <- data.frame(predicted=predict(model),
                   residual = residuals(model))
ggplot(tdat,aes(x=residual)) + geom_histogram(bins=20, color="black")
ggplot(tdat,aes(sample=residual)) + stat_qq() + stat_qq_line()


(model<-lme(tipping.points ~ 
               degree + nestedness  + a, random = ~1|web, data = new_df,
             na.action = "na.fail"))

summary(model)
anova(model)
dredge(model)




library(viridis)

new_df %>%
ggplot(aes(y=sd, x=degree, group=patches,color = factor(patches)))+
  geom_point(size=2.75, alpha=0.5)+
  stat_smooth(method = "loess", se=F)+
  theme_bw()+theme(legend.position = "bottom")+xlab("species degree")+
  ylab("AUC of species SD")+
  facet_wrap(.~factor(a))

new_df %>% 
  ggplot(aes(y=Autocorrelation, x=degree, group=patches,color = factor(patches)))+
  geom_point(size=2.75, alpha=0.5)+
  #stat_smooth(method = "loess", se=F)+
  theme_bw()+theme(legend.position = "bottom")+xlab("species degree")+
  ylab("AUC of species autocorrelation")+
  facet_wrap(.~factor(a))


new_df %>% 
  ggplot(aes(y=tipping.points, x=(degree), color = factor(patches)))+
  geom_point(size=2.75, alpha=0.5)+
 # stat_smooth(method = "loess", se=F,formula = y ~ x)+
  theme_bw()+theme(legend.position = "bottom")+xlab("species degree")+
  ylab("Mutualistic strength at collapse")+
  facet_wrap(.~factor(a))


#historgrams


(d1<-new_df %>% filter(degree < 12 ) %>% 
  ggplot(aes(x=tipping.points, fill = factor(patches)))+
  geom_histogram()+ggtitle(bquote("A) species degree:"< 12))+
  theme_cowplot()+xlab(" Tipping point")+
  facet_wrap(.~a))


(d3<-new_df %>% filter(degree > 12 ) %>% 
  ggplot(aes(x=tipping.points, fill = factor(patches)))+
  geom_histogram()+ggtitle(bquote("B) species degree" > 12))+
  theme_cowplot()+xlab(" Tipping point")+
  facet_wrap(.~a))


grid_arrange_shared_legend(d1,d3,nrow=2,ncol=1)

#averaging the tipping points
new_df_avg<-new_df %>% 
  group_by(web,nestedness, connectance, a, patches, No.species) %>% 
  summarise(mean_tipping_point = mean(tipping.points,na.rm=T),
            median_tipping_point = median(tipping.points),
            median_degree = median(degree),
            mean_degree = mean(degree),
            sd_tipping_point = sd(tipping.points,na.rm=T))


mod_1<-(lme(mean_tipping_point~ nestedness+ as.factor(a) + median_degree, random = ~1|web,
    data=new_df_avg))
summary(mod_1)



ggplot(new_df_avg,aes(x=connectance, y=median_tipping_point, color = factor(patches)))+
  geom_point()+
  theme_cowplot()+
  facet_wrap(.~factor(a))

