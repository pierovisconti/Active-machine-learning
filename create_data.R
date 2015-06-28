

####### CREATE DATASET #########
create_data <- function(xmax,ymax,distro,prevalence,meanx=xmax/2,meany=ymax/2,standx=0.1*xmax,standy=0.1*ymax,...){
# THIS READS
  # xmax,ymax: maximum coordinates
  # distro: the distribution of presence values: either normal or uniform
  # meanx and meany: the mean of the 2D gaussian coordinate distribution of the distribution
  # 
### DEFINE COORDINATES AND INITIALIZE OUTPUT
xcoord<-1:xmax
ycoord<-1:ymax


## RETURN SPECIES DISTRIBUTION
output<-switch(distro,
               uniform=uniform(size,prevalence,output=matrix(0,xmax,ymax)),
               normal=normal(xmax,ymax,prevalence,size,meanx,meany,standx,standy,output=matrix(0,xmax,ymax)),
               polynomial=polynomial(xmax,ymax,prevalence,...)
) 
return(output)
}

## DEFINE DISTRIBUTION FUNCTIONS
uniform <- function(size,prevalence,output=matrix(0,xmax,ymax)) {
  ind=sample(size,prevalence*(size))
  output[ind]=1
  return(output)
}
normal<- function(xmax,ymax,prevalence,size,meanx=xmax/2,meany=ymax/2,standx=0.1*xmax,standy=0.1*ymax,output=matrix(0,xmax,ymax)){
  px=as.matrix(dnorm(1:xmax,mean=meanx,sd=standx))
  py=as.matrix(dnorm(1:ymax,mean=meany,sd=standy))
  probs=apply(px,1,function(x) x*py)
  ind=sample(1:size,prevalence*size,prob=probs)
  output[ind]=1
  return(output)
}

polynomial<-function(xmax,ymax,prevalence,predictors,proc_error,interc,t1,t2,t3,p1,p2,p3){
  probs=1/(1+exp(-(interc+t1*predictors[,1]+t2*predictors[,1]^2+t3*predictors[,1]^3+
                     p1*predictors[,2]+p2*predictors[,2]^2+p3*predictors[,2]^3)))
output=matrix(0,xmax,ymax)
if (proc_error==TRUE){
ind=sample(1:size,round(prevalence*size),prob=probs^3)
} else {
  ind=order(probs,decreasing=TRUE)[1:round(prevalence*size)]
}
output[ind]=1
return(output)
}


