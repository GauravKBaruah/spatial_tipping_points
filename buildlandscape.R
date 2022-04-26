

##############################################################
### PatchValue function
### Function takes two  parameters:
###     N: number of patches 
###     stdvalues: standard dev of patch values
###     valid parameters --> N integer > 0, stdvalues > 0 
### The function returns a vector with random independent patch values. The average patch value is fixed to one and the standar deviation is equal to stdvalues. Patch values are distributed as a beta distribution (rescaled to have support between 0 and 2)
PatchValue <- function( N , stdvalues ){
  # the distribution of patch values is a (rescaled) beta distribution. Change the following line to change the distribution
  if ( stdvalues == 0. ){
    PatchValues <- rep( 1. , N )
  } else{
    PatchValues <- 2. * rbeta(N, 0.5*(1./stdvalues^2 - 1), .5*(1./stdvalues^2 -1) ) 
  }
  return( PatchValues )
}





##############################################################
### distperiodic function
### Function takes two parameters:
###     x: vector of coordinates: N patches and d coordinates (d = number of dimensions) 
###     L: the landscape has side L (and volume L^d)
###     valid parameters --> components of x > 0 and < L, L > 0
### the function returns a matrix with the distances between all the patches. The distance is computed using periodic boundary conditions 
distperiodic <- function(x, L){
  N <- dim(x)[1] # number of patches
  D <- dim(x)[2] # number of spatial dimensions
  Distances <- matrix(0, N, N) # initialize the matrix of distances
  for ( dd in 1:D ){
    coordinate <- x[,dd]
    twopoints <- expand.grid(coordinate, coordinate) # Create a data frame from all combinations of coordinates
    contribution <- matrix(abs(twopoints[,1] - twopoints[,2]), N , N )
    contribution[contribution > (0.5 * L)] <- L - contribution[contribution > (0.5 * L)]
    Distances <- Distances + contribution^2 # add the contribution od distance fiven by componend dd
  }
  return(sqrt(Distances))
}


##############################################################
### BuildRndLandscape function
### Function takes two parameters:
###     N: number of patches 
###     d: number of dimension
###     w: dispersal kernel
###     L: the landscape has side L (and volume L^d) [default = 1]
###     Kernel: name of a kernel in Kernels.R
###     sigmaA: std of patch values distribution
###     Bounds: can be open [default] or periodic 
###     valid parameters --> N positive integer, d positive integer, L > 0, xi > 0, FUN a function declared in Kernels.R, sigmaA > 0, Bounds = 'Open' or 'Periodic' 
### the function returns the community matrix M
BuildRndLandscape <- function( N, d, L = 1, w, Coordinates, FUN = Exponential, sigmaA , Bounds = 'Open'){
  
  
  ## Calculate the distances between patches
  if ( Bounds == "Periodic" ){
    Distances <- as.matrix(distperiodic(Coordinates,L))
  }
  else{
    Distances <- as.matrix(dist(Coordinates))
  }
  
  
  Asizes <- PatchValue(N, sigmaA ) # vector of patch values
  
  SizeMat <- Asizes %*% t( Asizes )
  
  M <- FUN( Distances, w ) * SizeMat
  
  return( M )
  
}


##############################################################
### plot landscape function
### Function takes two parameters:
###     N: number of patches 
###     d: number of dimension
###     w: dispersal kernel
###     L: the landscape has side L (and volume L^d) [default = 1]
###     Kernel: name of a kernel in Kernels.R
###     sigmaA: std of patch values distribution
###     Bounds: can be open [default] or periodic 
###     valid parameters --> N positive integer, d positive integer, L > 0, xi > 0, FUN a function declared in Kernels.R, sigmaA > 0, Bounds = 'Open' or 'Periodic' 
### the function returns the community matrix M

plot_net<-function(mat,Coordinates){ 
  net<-graph_from_adjacency_matrix(mat, mode="undirected",
                                   weighted = TRUE, diag=TRUE)
  colrs<- c("gray50", "tomato", "firebrick")
  V(net)$colr <- colrs[V(net)$media.type]
  deg<-degree(net,mode="all")
  V(net)$size <- 7  #deg*0.5
  V(net)$label <- NA
  E(net)$width <- E(net)$weight*6.5
  E(net)$arrow.size <-2
  E(net)$edge.color < "firebrick"
  plot(net, layout = as.matrix(Coordinates), edge.color = "tomato")
  
}

## build fragmented dispersal landscape
#patches_removed= vector of patches removed
# d = dimensions which is 2
# L = length of the x-y laandscape : which is 1 here
# w = dispersal length
#Coordinates = coordinates in a x-y landscape
# FUN = dispersal funciton used

frag_disp<-function(patches_removed, L = 1, w, Coordinates, FUN = Exponential, sigmaA , Bounds = 'Open'){
  ## Calculate the distances between patches
  
  if ( patches_removed == 0) { newCoordinates <- Coordinates
  } else{ patches_rm_id<-sample(1:no.of.patches, patches_removed, replace = F)
    newCoordinates<-Coordinates[-patches_rm_id,]
    }
  
  N <- nrow(newCoordinates) # no. of new patches remaining
  
  if ( Bounds == "Periodic" ){
    Distances <- as.matrix(distperiodic(newCoordinates,L))
  }
  else{
    Distances <- as.matrix(dist(newCoordinates))
  }
  
  
  Asizes <- PatchValue(N, sigmaA ) # vector of patch values
  
  SizeMat <- Asizes %*% t( Asizes )
  
  M <- FUN( Distances, w ) * SizeMat
  
  return( M )
  
  
}
