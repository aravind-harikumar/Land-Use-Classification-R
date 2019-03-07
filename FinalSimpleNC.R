rm(list=ls(all=TRUE))
require(rgdal)
require(mvtnorm)
require(rgl)
require(scatterplot3d)
require(R.utils)

# Read the input image and extrat the dimensionality details
Path <- "F:\\Aravind\\RProgram\\Data"
str_name<-"abc.tif"
imageObj <- readGDAL(paste(Path , "\\", str_name, sep=""))
imageObjArray <- as.array(imageObj)
d <- dim(imageObjArray)
Bands <-d[3]
N <- d[2]
M <- d[1]
n <- M*N

#Fuzzification Factor
m =2.2

#Delta value
Delta <- 10000

# number of classes
Ncl <- 6

# Maximum DN Value
maxDN <- 255

# Randomly initialize membership values
MembValArray = array(0,c(M,N,Ncl+1))

# Assign mean values for each class and initialization of cluster centers
MeanClassVal <- array(0,c(Bands,Ncl))
MeanClassVal[,] <- c(c(83,39.77,102.45,55.4),c(76.27,36,103.83,54.77), c(77.78,39.15,84.15,47.31), c(113.22,82.27,101.05,125.05),c(91.88,54.88,56,71.88), c(70.76,34,27.03,16))

# Function to find initial optimal memebrship values for all classes and pixels
getAllUValues <- function(){
  Md <- array(0,c(Ncl,n))
  Mdk <- array(0,c(Ncl,Ncl,n))
  Mdy <- array(0,c(Ncl,1,n))
  
  aa = array(imageObjArray, c(n,Bands))
  bb = t(MeanClassVal)
  
  for(k in 1:Ncl)
    Md[k,] <- sqrt(rowSums((t(t(aa)-bb[k,]))^2,1))
  
  for(l in 1:Ncl)
    for (k in 1:Ncl){
      Mdk[k,l,] <- (Md[l,]/Md[k,])^(2/(m-1)) + (Md[l,]/Delta)^(2/(m-1))      
    }
  
  for(l in 1:1)
    for (k in 1:Ncl){
        Mdy[k,l,] <- (Delta/Md[k,])^(2/(m-1))
    }
  
  NoiseClsVector <- t(1.0/(colSums(Mdy,1)+1))
  ClassVectors <- t(1.0/colSums(Mdk,1))
  
  m1=array(0,c(M,N,Ncl+1))
  
  m1[,,1] <- matrix(ClassVectors[,1], nrow=M, ncol=N)
  m1[,,2] <- matrix(ClassVectors[,2], nrow=M, ncol=N)
  m1[,,3] <- matrix(ClassVectors[,3], nrow=M, ncol=N)
  m1[,,4] <- matrix(ClassVectors[,4], nrow=M, ncol=N)
  m1[,,5] <- matrix(ClassVectors[,5], nrow=M, ncol=N)
  m1[,,6] <- matrix(ClassVectors[,6], nrow=M, ncol=N)
  m1[,,7] <- matrix(NoiseClsVector, nrow=M, ncol=N)
  
  return(m1)  
} 

# # Function to find initial optimal cluster centers for a specific class
Vvalues <- function (cls){
  
  MemValueVector <- array(rep(MembValArray[,,cls]^m,times=Bands) , dim=c(M,N,Bands))
  cc <- MemValueVector * imageObjArray
  
  aaa <- c(sum(cc[,,1]), sum(cc[,,2]), sum(cc[,,3]), sum(cc[,,4]))
  return(aaa/(sum(MemValueVector)/4))
}

# Function to find initial optimal cluster centers for all classes 
getAllVvalues <- function (){
  for (cl in 1: Ncl)
    MeanClassVal[,cl] <- Vvalues(cl)
  return(MeanClassVal)
}

# Function to stretch the histogram for better results 
hist_stretch<-function(data)
{
  cur.lim<-quantile(data,c(0.025,0.975),na.rm=TRUE)
  data<-pmax(cur.lim[1],pmin(cur.lim[2],data))
  data<-floor(255*(data-cur.lim[1])/(cur.lim[2]-cur.lim[1]))
  data[is.na(data)]<-0
  data
}

# Objective function of NC
NC <- function(){
  firstTerm <- 0
  secondTerm <- 0
  
  for (i in 1:M)
    for (j in 1:N){
      for(cl in 1:Ncl) {    
        firstTerm <- firstTerm + ((MembValArray[i,j,cl])^m)*(sqrt(sum((imageObjArray[i,j,] - MeanClassVal[,cl])^2)))       
      }
      secondTerm <- secondTerm + ((MembValArray[i,j,Ncl+1])^m)*Delta
    }
  return(firstTerm+secondTerm)
}

temp <- MembValArray

objfunvalold <- 0
objfunvalnew <- 0
Objerr <- 0

# Code to create a SpatialGridDataFrame
gsd <- 56
xyoffset <- c(271761.75, 5822778.0)
xyoffset <- xyoffset - 0.5*c(gsd,gsd)
Refdata <-data.frame(as.vector(imageObjArray[,,1]))
names(Refdata) <- "class"
refgrid <- GridTopology(cellcentre.offset=xyoffset,cellsize=c(gsd,gsd),cells.dim=c(M,N))
Ref <- SpatialGridDataFrame(refgrid, Refdata, proj4string = "+proj=utm +zone=44 +datum=WGS84 +units=m +no_defs +ellps=WGS84
                            +towgs84=0,0,0")
Ref$class <-  hist_stretch(Ref$class)
windows()

#The beloe code implements the main process of NC algorithm

iter<-0
for (iter in 1:1){
  ptm <- proc.time()
  
  #get U values (i.e membership values)
  MembValArray <- getAllUValues()
  MembValArray[MembValArray==Inf] <-0
    
  objfunvalnew <- NC() 

  # Divergence/Convergence value of objective function (i.e current iteration value - previous iteration value)
  Objerr <- max(abs(objfunvalold-objfunvalnew))
  objfunvalold <- objfunvalnew
  print(Objerr)
  
  # Loop exit condition
  if(Objerr<50)
    break  
  
  MLC <- Ref
  MLC$class <- array(0,n)
  
  # Code for display
  par(mfrow=c(2,4))  
  for(i in 1:(Ncl+1)){
    F <- as.vector(MembValArray[,,i])
    MLC$class <- round((F/max(F)) * 255) 
    image(MLC, col=gray((0:255)/255), axes=TRUE)
    title(paste({"Class "}, {i-1} ,{" - Iteration "},iter, sep=""))
  }
  
  #get V values (i.e cluster centers/ mean vector for each class)
  MeanClassVal <- getAllVvalues()
  
  print(proc.time() - ptm)
}