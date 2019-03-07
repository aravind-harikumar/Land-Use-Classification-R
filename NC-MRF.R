rm(list=ls(all=TRUE))
require(rgdal)
require(mvtnorm)
require(rgl)
require(scatterplot3d)
require(R.utils)

# Read the input image & extrat the dimensionality details
Path <- "F:\\Aravind\\RProgram\\Data"
str_name<-"abc.tif"
imageObj <- readGDAL(paste(Path , "\\", str_name, sep=""))
imageObjArray <- as.array(imageObj)
d <- dim(imageObjArray)
Bands <-d[3]
N <- d[2]
M <- d[1]
n <- M*N

# Actual window size is 2*WSize+1
WSize  <- 1

# Max number of neighbours
Nn  <- (WSize*2+1)^2-1

#Fuzzification Factor
m =2.2

#Delta value
Delta <- 1000

# number of classes
Ncl <- 6

# The weightage to spatial and spectral components
lambda <- 0.9

# Gamma value
gamma <- 0.5

# Beeta value
beeta <- 2

#Selected Prior
prior = "DA2"  #  Options: "SA", "DA1", "DA2", "DA3", "DA4"

# Maximum DN Value
maxDN <- 255

Neigh_Coord  <- array(0, c(M, N, Bands))
Weight  	<- array(0, c(2*WSize+1, 2*WSize+1))

# Randomly initialize membership values
MembValArray = array(runif(N*M*(Ncl+1)),c(M,N,(Ncl+1)))

# Assign mean values for each class (Since it is supervised method)
MeanClassVal <- array(0,c(Bands,Ncl))
MeanClassVal[,] <- c(c(83,39.77,102.45,55.4),c(76.27,36,103.83,54.77), c(77.78,39.15,84.15,47.31), c(113.22,82.27,101.05,125.05),c(91.88,54.88,56,71.88), c(70.76,34,27.03,16))

# Function to find initial optimal memebrship values for all classes and pixels
getInitMembValues <- function(){
  Md <- array(0,c(Ncl,n))
  Mdk <- array(0,c(Ncl,Ncl,n))
  Mdy <- array(0,c(Ncl,1,n))
  
  # Vectorized bands for easy computaion
  vectorizedBands = array(imageObjArray, c(n,Bands))
  meanVector = t(MeanClassVal)
  
  # Computaion steps
  for(k in 1:Ncl)
    Md[k,] <- sqrt(rowSums((t(t(vectorizedBands)-meanVector[k,]))^2,1))
  
  for(l in 1:Ncl)
    for (k in 1:Ncl){
      Mdk[k,l,] <- (Md[l,]/Md[k,])^(2/(m-1)) + (Md[l,]/Delta)^(2/(m-1))      
    }
  
  for(l in 1:1)
    for (k in 1:Ncl){
      Mdy[k,l,] <- (Delta/Md[k,])^(2/(m-1))
    }
  
  # New membership values in Vector form 
  NoiseClassVector <- t(1.0/(colSums(Mdy,1)+1))
  ClassVectors <- t(1.0/colSums(Mdk,1))
  
  #Convert back to Matrix of dimention  M*N*Bands (Rows/Column is set to M/N to match with 'MembValArray')
  m1=array(0,c(M,N,Ncl+1))
  for (bnd in 1:Ncl){
    m1[,,bnd] <- matrix(ClassVectors[,bnd], nrow=M, ncol=N)
  }  
  m1[,,(Ncl+1)] <- matrix(NoiseClassVector, nrow=M, ncol=N) # For new Noise Membership (Equation is different)
  
  return(m1)  
}

# Get Updated Membership values for each class.
UNew <- function(UMat){
  
  tempBands = array(UMat, c(n,Ncl+1)) 
  tempBands1 = array(UMat, c(n,Ncl+1)) 
  
  # Code for generating membership value from shuffling membership values (not within bands but across bands) of input UMat 
  for (J in 1:Ncl+1) 
    for (val in 1:n){ 
      tempBands1[val,J]= tempBands[val,(round(runif(1,min=1,max=977))%%J)+1]
    }
  
  #Code for generating random membership value generation
  
#     for (val in 1:n){ 
#       tempBands1[val,]= sample_frac()  
#     }
  
  m1=array(0,c(M,N,Ncl+1))
  
  for(bnd in 1:(Ncl+1))
    m1[,,bnd] <- matrix(tempBands1[,bnd], nrow=M, ncol=N)
  
  return(m1)   
}

# Find the L2 norm of a vector
L2distObjFunc <- function(rowNo, colNo, classNo){  
  dist <-sqrt(sum((imageObjArray[rowNo,colNo,] - MeanClassVal[,classNo])^2))
  return(dist)
}

# Noise Clustering Objective Function
NC <- function(i,j,cl){
  totVal <- 0  
  totVal <- totVal + ((f[i,j,cl])^m)*L2distObjFunc(i,j,cl) + ((f[i,j,(Ncl+1)])^m)*Delta
  return(totVal)
}

# Function assigning weights in the neighbourhood
Fw <- function(a,b){  
  val <- a^2 + b^2
  val <- 1 / val
  val <- val^(0.5)
  val[val==Inf]<-0
  return(val)
}

sample_frac<-function()
{
  val <- array(0,Ncl+1)  
  k <- round(runif(1,min=0,max=(Ncl+1))+0.5)  
  val[k] <- runif(1,min=1/(Ncl+1),max=1.0)  
  k_rest <- (1:(Ncl+1))[-k]  
  f_lim <- 1-val[k]  
  f_rest <- runif(Ncl,min=0,max=f_lim)
  s<-sum(f_rest)  
  while(s==0)
  {
    f_rest <- runif(Ncl,min=0,max=1)
    s<-sum(f_rest)
  }  
  val[k_rest] <- f_rest*f_lim/s  
  return(val)
}

# for(k in 1:(2*WSize+1))
#   for(l in 1:(2*WSize+1))
#   {
#     Weight[k, l] <- Fw(k-(WSize+1),l-(WSize+1))
#   }
# 
# Weight <- Weight/ sum(Weight)

for(i in 1:M)
  for(j in 1:N)
  {
    imin <- i - WSize
    imax <- i + WSize
    jmin <- j - WSize
    jmax <- j + WSize
    
    if(imin<1)  imin <-1
    if(imax>M) imax <-M
    if(jmin<1)  jmin <-1
    if(jmax>N) jmax <-N
    
    Neigh_Coord[i, j, ] <- c(imin,imax,jmin,jmax)
  }

Uprior <- function(i,j,cl){
  val <- 0
  
  f1 <- f[((Neigh_Coord[i,j,1]):(Neigh_Coord[i,j,2])),((Neigh_Coord[i,j,3]):(Neigh_Coord[i,j,4])),cl]
  #W1 <- Weight[(Neigh_Coord[i,j,1]-i+1+WSize):(Neigh_Coord[i,j,2]-i+1+WSize),(Neigh_Coord[i,j,3]-j+1+WSize):(Neigh_Coord[i,j,4]-j+1+WSize)]
  
  f0 <- (f1 - f[i,j,cl])^2
  
  for (ct in 1:(dim(f0)[1]*dim(f0)[2])){
  #  if(ct == (round((dim(f0)[1]*dim(f0)[2])/2)+1))
   #   next
    val <- val + (beeta * (f0[ct]))
  }

  return(val)
}

DA1prior <- function(i,j,cl){
  val <- 0
  
  f1 <- f[((Neigh_Coord[i,j,1]):(Neigh_Coord[i,j,2])),((Neigh_Coord[i,j,3]):(Neigh_Coord[i,j,4])),cl]
  #W1 <- Weight[(Neigh_Coord[i,j,1]-i+1+WSize):(Neigh_Coord[i,j,2]-i+1+WSize),(Neigh_Coord[i,j,3]-j+1+WSize):(Neigh_Coord[i,j,4]-j+1+WSize)]
  
  f0 <- (f1 - f[i,j,cl])^2
  
  for (ct in 1:(dim(f0)[1]*dim(f0)[2])){
    val <- val + (-gamma * (exp((-(f0[ct])/gamma))))  
  }
  return(val)
}

DA2prior <- function(i,j,cl){
  val <- 0
  
  f1 <- f[((Neigh_Coord[i,j,1]):(Neigh_Coord[i,j,2])),((Neigh_Coord[i,j,3]):(Neigh_Coord[i,j,4])),cl]
  #W1 <- Weight[(Neigh_Coord[i,j,1]-i+1+WSize):(Neigh_Coord[i,j,2]-i+1+WSize),(Neigh_Coord[i,j,3]-j+1+WSize):(Neigh_Coord[i,j,4]-j+1+WSize)]
  
  f0 <- (f1 - f[i,j,cl])^2
  
 for (ct in 1:(dim(f0)[1]*dim(f0)[2])){
#    if(ct == (round((dim(f0)[1]*dim(f0)[2])/2)+1))
#      next
    val <- val + (-gamma/(1 + ((f0[ct])/gamma)))
  }
  return(val)
}

DA3prior <- function(i,j,cl){
  val <- 0
  
  f1 <- f[((Neigh_Coord[i,j,1]):(Neigh_Coord[i,j,2])),((Neigh_Coord[i,j,3]):(Neigh_Coord[i,j,4])),cl]
  #W1 <- Weight[(Neigh_Coord[i,j,1]-i+1+WSize):(Neigh_Coord[i,j,2]-i+1+WSize),(Neigh_Coord[i,j,3]-j+1+WSize):(Neigh_Coord[i,j,4]-j+1+WSize)]
  
  f0 <- (f1 - f[i,j,cl])^2
  
  for (ct in 1:(dim(f0)[1]*dim(f0)[2])){
    val <- val + (gamma * log(1 + ((f0[ct])/gamma)))
  }
  
  return(val)
}

DA4prior <- function(i,j,cl){
  val <- 0
  
  f1 <- f[((Neigh_Coord[i,j,1]):(Neigh_Coord[i,j,2])),((Neigh_Coord[i,j,3]):(Neigh_Coord[i,j,4])),cl]
  #W1 <- Weight[(Neigh_Coord[i,j,1]-i+1+WSize):(Neigh_Coord[i,j,2]-i+1+WSize),(Neigh_Coord[i,j,3]-j+1+WSize):(Neigh_Coord[i,j,4]-j+1+WSize)]
  
  f0<- (f1 - f[i,j,cl])
  
  for (ct in 1:(dim(f0)[1]*dim(f0)[2])){
    val <- val + (gamma*abs((f0[ct])) - (gamma^2)*log(1 + (abs((f0[ct]))/gamma)))
  }
  
  return(val)
}

U <- function(i,j,clv){  
  val <- 0  
  if(prior == "SA")
    val <-  (1.0-lambda) * NC(i,j,clv) + lambda * Uprior(i,j,clv)
  else if(prior == "DA1")
    val <-  (1.0-lambda) * NC(i,j,clv) + lambda * DA1prior(i,j,clv)
  else if(prior == "DA2")
    val <-  (1.0-lambda) * NC(i,j,clv) + lambda * DA2prior(i,j,clv)
  else if(prior == "DA3")
    val <-  (1.0-lambda) * NC(i,j,clv) + lambda * DA3prior(i,j,clv)
  else
    val <-  (1.0-lambda) * NC(i,j,clv) + lambda * DA4prior(i,j,clv)

  return(val)
}

#================================================================================================
# Block 5:  Energy optimisation with simulated annealing
#================================================================================================

# Maximum allowed number of iterations
Niter <- 20

# SA cooling schedule parameters
T0 <- 2.0

# Convergence criterion for SA
min_acc_thr <- 0.1*10^(-2)*n*Ncl

T <- T0

par(mfrow=c(1,1))
stop_crit <- 0

# Initialize SpatialGridDataFrame
gsd <- 20
xyoffset <- c(271761.75, 5822778.0)
xyoffset <- xyoffset - 0.5*c(gsd,gsd)
#Refdata <- data.frame(as.vector(array(0,c(M,N))))
Refdata <-data.frame(as.vector(imageObjArray[,,1]))
names(Refdata) <- "MembVal"
refgrid <- GridTopology(cellcentre.offset=xyoffset,cellsize=c(gsd,gsd),cells.dim=c(M,N))
Ref <- SpatialGridDataFrame(refgrid, Refdata, proj4string = "+proj=utm +zone=44 +datum=WGS84 +units=m +no_defs +ellps=WGS84
                            +towgs84=0,0,0")

windows()

# Initialize the NC Array
NCA <- Ref
NCA$MembVal <- array(0,n)

# Get the inital memebership values
f <- getInitMembValues()
#f[f==Inf] <-0

MembValArray <- f
UpdatedUVal=UNew(MembValArray)

for(iter in 1:Niter)
{
  upd_count <- 0
  
  for(cl in 1:Ncl)
  {
    if(iter == 1)
      next
    
    for(i in 1:M)
    {
      for(j in 1:N)
      {        
        f_update <-  UpdatedUVal[i,j,cl]
        ft <- f[i,j,cl]
        
        if(f_update!=ft)
        {
          u1 <- U(i,j,cl)
          f[i,j,cl] <- f_update
          
          u2 <- U(i,j,cl)
          
          u1 <- u2-u1
          
          if(T!=0)
          {    
            if(u1>0)
            {              
              f[i,j,cl] <- ft
            }
            else
            {
              u1 <- exp(-u1/T)
              xi <- runif(1, min=0, max=1)
              
              if(xi>u1)
              {              
                f[i,j,cl] <- ft
              }
              else 
                upd_count<-upd_count+1
            }
          }              
        }   
      }
    }
  }
  
  if(iter!= 1 && upd_count<=min_acc_thr) # Iter 1 is just output of Simple NC i.e Initial Membership values
    break
  
  if(iter!= 1) # Iter 1 is just output of Simple NC i.e Initial Membership values
    print(upd_count)
  
  # Update Temperature
  T <- (log(1+iter)/log(2+iter))*T0
  
  # Print the NC-MRF output for each class. (i.e Membership values for each class)
  par(mfrow=c(2,4))
  
  for(i in 1:(Ncl+1)){
    F <- as.vector(f[,,i])
    NCA$MembVal <- round((F/max(F)) * 255) 
    image(NCA, col=gray((0:255)/255), axes=TRUE)
    
    if(iter==1)
      title(paste({"Class "}, {i-1}, {" Init Membership Values"}))
    else
      title(paste({"Class "}, {i-1}, {" - Iteration "},(iter-1), sep=""))  # First iter is just skipped over
  }
  
}