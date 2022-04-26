
##############################################################
### Dispersal kernels functional forms
### Function takes two (or more) parameters:
###     distance: distance between two patches 
###     xi: dispersal length
###     alpha: shape of the beta function (in the case of beta kernels)
###     beta: shape of the beta function (in the case of beta kernels)
###     valid parameters --> distance > 0, xi > 0, alpha > 0, beta > 0
### distance can be either a number, a list or a matrix. 
### The function returns either a number, a list or a matrix (depending on the type of distance),
### whose value(s) is the dispersal rate at a given distance

thompson_etal<-function(distance, w){
  
  M<-exp(-distance/w)
  M<-M-diag(nrow(M))
  M[M==1]<-0
  disp_matrix<-apply(M, 1, function(x) x/sum(x))
  
  return(disp_matrix)
}


