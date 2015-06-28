rm(list=ls())
setwd('C:/Users/a-pierov/Dropbox/Visconti_Joppa/bestpoint/scripts')
######################### DEFINE LANDSCAPE ######################
xmax=50
ymax=50
size=xmax*ymax
data=matrix(NA,xmax,ymax) # prior data
numpred=2
samples=as.data.frame(matrix(nrow=0,ncol=3+numpred)) ## initial list of samples
namespred=c('temp','prec')
names(samples)=c("time","location","PA",namespred )
coords=1:(xmax*ymax) ## locations available for sampling

## CREATE INPUT VARIABLES FOR FULL LANDSCAPE MODEL
xcoord=as.data.frame(matrix(1:xmax,nrow=xmax*ymax,byrow=TRUE))
ycoord=as.data.frame(c(matrix(1:ymax,nrow=ymax,ncol=xmax,byrow=TRUE)))
predictors=as.data.frame(cbind(xcoord,ycoord))
#predictors[180:225,]=10
names(predictors)=namespred

## CREATE LANDSCAPE
source('create_data.R')
source('sampling_methods.R')
prevalence=0.2
#interc=0.5;t1=0.2;t2=-0.05;t3=0;p1=0.015;p2=0.023;p3=0;
interc=0;t1=0.05;t2=-0.0025;t3=0;p1=0.030;p2=-0.002;p3=0;
landscape<-create_data(xmax,ymax,'normal',prevalence,7,7,0.5,0.5)# create landscape 
landscape<-create_data(xmax,ymax,'polynomial',prevalence,predictors=predictors,proc_error=FALSE,interc=interc,t1=t1,t2=t2,t3=t3,p1=p1,p2=p2,p3=p3)
landscape
labels=factor(matrix(c(landscape),nrow=xmax*ymax))


## DEFINE SAMPLING EFFORT
nsamples=20
npoints=1 # HOW MANY POINTS FOR EACH TIME STEP

## "pre-allocate" an empty list of length nsamples and AUC values
models=vector("list",nsamples) 
AUCvalues=matrix(,1,nsamples)
TSSvalues=matrix(,1,nsamples)
threshold=matrix(,1,nsamples)
prior_data=list(samples=samples,data=data,coords=coords)
results=list(prior_data,models=models,AUC=AUCvalues,TSS=TSSvalues,threshold=threshold) ## CHECK IF THE THRESHOLD WORKS

clus=clara(predictors,9,sampsize=225)
clusters=clus$clustering
clusplot(clus)
palette(rainbow(9))
par(mfrow=c(1,1))
plot(predictors[,1],predictors[,2],type='p',pch=22, col=clusters,bg=clus$clustering,cex=9)


library(AUC)
library(GRaF)
library(cluster)



#### ################################# TEST    ###############################
## GENERATE PRIOR DATA
formula='PA ~ temp +prec +I(temp^2) + I(prec^2) + I(temp^3) + I(prec^3)'

#formula='PA ~ -1 temp +prec'
prior_data=list(samples=samples,data=data,coords=coords)
prior_data=sampling('strat',time_stamp=1,prior_data,predictors,'GRAF',formula,labels,landscape,npoints=10,predictions=0,threshold=0.5)
results_prior=active_learning(nsamples=2,method='uniform',modeltype='GRAF',formula,labels,prior_data,predictors,landscape,npoints=5)
prior_data=results_prior[[1]]
## TRY  RANDOM SAMPLING
nreps=10
nsamples=20
AUCunif_graf=matrix(nrow=nsamples,ncol=nreps)
TSSunif_graf=matrix(nrow=nsamples,ncol=nreps)
AUCstrat_graf=matrix(nrow=nsamples,ncol=nreps)
TSSstrat_graf=matrix(nrow=nsamples,ncol=nreps)
AUCunc_graf=matrix(nrow=nsamples,ncol=nreps)
TSSunc_graf=matrix(nrow=nsamples,ncol=nreps)
AUCma_graf=matrix(nrow=nsamples,ncol=nreps)
TSSma_graf=matrix(nrow=nsamples,ncol=nreps)
for (i in 1:nreps) {
  
  print (paste('------------------------',i,'---------------------'))  
  
  indicator=0
  while (indicator<2 | indicator==6) {
    prior_data=list(samples=samples,data=data,coords=coords)
    prior_data=sampling('uniform',time_stamp=1,prior_data,predictors,'GRAF',formula,labels,landscape,npoints=6,predictions,threshold=0.5)
    indicator=sum(prior_data$samples[,3])
    #results_uniform=tryCatch(active_learning(nsamples=10,method='uniform',modeltype='GRAF',formula,labels=labels,prior_data=prior_data,predictors=predictors,landscape=landscape,
                                         #    npoints=1),error=function(c){ message("----------Try again!--------"); return('yes')})
  }

  results_unif=active_learning(nsamples=nsamples,method='uniform',modeltype='GRAF',formula,labels=labels,prior_data=prior_data,predictors=predictors,landscape=landscape,
                                npoints=1)
  AUCunif_graf[,i]=results_unif$AUC
  TSSunif_graf[,i]=results_unif$TSS[1:nsamples]
#   
  
  results_strat=active_learning(nsamples=nsamples,method='strat',modeltype='GRAF',formula,labels=labels,prior_data=prior_data,predictors=predictors,landscape=landscape,
                                 npoints=1)
   #results_unc_pred=active_learning(nsamples=15,'uncert_pred',modeltype='GRAF',formula,labels=labels,prior_data=prior_data,predictors=predictors,landscape=landscape,
   #                            npoints=1)
AUCstrat_graf[,i]=results_strat$AUC
TSSstrat_graf[,i]=results_strat$TSS[1:nsamples]

#     AUCunc_graf[,i]=results_unc_pred$AUC
#     TSSunc_graf[,i]=results_unc_pred$TSS[1:15]
#   
 results_maxent=active_learning(nsamples=nsamples,'maxent',modeltype='GRAF',formula,labels=labels,prior_data=prior_data,predictors=predictors,landscape=landscape,
                                  npoints=1)
#   #results_MIMC=active_learning(nsamples=10,'MIMC',modeltype='GRAF',formula,labels=labels,prior_data=prior_data,predictors=predictors,landscape=landscape,
#   #    npoints=1)
   AUCma_graf[,i]=results_maxent$AUC
   TSSma_graf[,i]=results_maxent$TSS[1:nsamples]
}

# UNIFORM SAMPLING
rowMeans(AUCunif_graf[,1:nreps])
rowMeans(TSSunif_graf[,1:nreps])


# UNCERTAINTY SAMPLING
rowMeans(AUCunc_graf[,1:nreps])
rowMeans(TSSunc_graf[,1:nreps])

## TRY MAXENT SAMPLING (FROM THE FAMILY OF UNCERTAINTY SAMPLING METHODS)
rowMeans(AUCma_graf[,1:nreps])
rowMeans(TSSma_graf[,1:nreps])

# STRATIFIED SAMPLING
rowMeans(AUCstrat_graf[,1:nreps])
rowMeans(TSSstrat_graf[,1:nreps])



#results_MIMC=active_learning(nsamples=2,'MIMC',modeltype='GRAF',formula,labels,prior_data,predictors,landscape,npoints=1)
#results_MIMC$AUC
#results_MIMC$TSS



# locs=prior_data_uniform$samples[1:5,2]
# data_maxent=data
# data_maxent[locs]=landscape[locs]
# threshold=results_uniform$threshold[1]
# models=results_uniform[[2]]
# input_pred=as.data.frame(predict(models[[1]],newdata=predictors,type='response'))
# prior_data_maxent=list(samples=prior_data_uniform$samples[1:5,],data=data_maxent,coords=coords[-locs])
# results_maxent=active_learning(nsamples=4,'maxent',prior_data_maxent,init=1,landscape,npoints=5,threshold,input_pred)
## FIT MODEL
# model_poly=glm(PA~-1+poly(xcoord,2)+poly(ycoord,2), family=binomial,data=dataset)
# summary(model_poly)






