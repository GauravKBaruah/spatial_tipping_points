
#code borrowed from Grilli et al 2015, Plos computational biology, 

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


    Distances <- as.matrix(dist(Coordinates))
  
  
  
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
