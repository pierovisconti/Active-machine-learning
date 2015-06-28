## DEFINE SAMPLING METHODS

uniform_samp <- function(xmax,ymax,coords,npoints,data=matrix(0,xmax,ymax),range,...) {
  ind=sample(coords,npoints)
  data[ind]=range[ind]
  return(data)
}  

stratified_samp <- function(coords,npoints,clusters,data,landscape) {
  for (i in 1:npoints){
  available<-clusters[coords]
  num_obs<-tapply(available,available,length)
  ind1<-which(num_obs==max(num_obs))
  if (length(ind1)==1) {
    ind2=ind1}else {ind2<-sample(c(ind1),1)}
  ind3<-coords[which(available==ind2)]
  ind<-sample(ind3,1)
  data[ind]<-landscape[ind]
  coords<-coords[-which(coords==ind)]
  }
  return(data)
}
  
  
  
 

# IT CHOOSE THE POINT WITH HIGHEST PREDICTED  ENTROPY: P*LOG(P), I.E. HIGHEST UNCERTAINTY
maxent_samp<-function(coords,prevalence,predictions,threshold,npoints,data=matrix(0,xmax,ymax),range){
  # THIS SELECTS THE POINT WHOSE MEDIAN POSTERIOR PREDICTION IS CLOSEST TO THE THRESHOLD FOR P/A
  dat=cbind(coords,predictions[coords,1])
  #npoints.pos=round(npoints/2)
  #npoints.neg=npoints-npoints.pos
    dist.threshold=abs(dat[,2]-threshold)
#   
#   # WHEN YOU ONLY CARE ABOUT UNCERTAINTY AND NOT TO HAVE A TRAINING DATASET OF THE RIGHT PREVALENCE USE THIS.
  ordered.dt=order(abs(dist.threshold))
  ind=coords[ordered.dt[1:npoints]]
  
  #  # THE FOLLOWING TRIES TO HAVE A TRAINING DATASET OF THE RIGHT PREVALENCE. 
#  
#   if (threshold>prevalence) {
#     dist.threshold.neg=dist.threshold; dist.threshold.neg[which(dist.threshold>=0)]=10
#     ordered.dt.neg=order(abs(dist.threshold.neg))
#     ind=coords[ordered.dt.neg[1:npoints]]
#   } else {
#     dist.threshold.pos=dist.threshold; dist.threshold.pos[which(dist.threshold<0)]=10
#     ordered.dt.pos=order(abs(dist.threshold.pos))
#     ind=coords[ordered.dt.pos[1:npoints]]
#   }
#    

  data[ind]=range[ind]
  return(data)
}

uncert_samp<-function(coords,predictions,npoints,data=matrix(0,xmax,ymax),range){
  # THIS SELECTS THE POINT WITH THE WIDEST CI IN PREDICTION
  unc=predictions[coords,3]-predictions[coords,2]
  ordered=order(unc,decreasing=TRUE)
  ind=coords[ordered[1:npoints]]
  data[ind]=range[ind]
  return(data)
}

# THIS CALCULATE P*LOG(P) OF ALL UNLABELLED INSTANCES. NEEDED FOR MCMI ALGORITHM. NEED TO INCORPORATE CT here
calc_entropy_unlabelled <- function (dataset,formula,alldata,i,modeltype,labels){
  dat=rbind(dataset,as.matrix(alldata[i,]))
  model<-switch(modeltype,
                GRAF=graf(dat$PA,dat[,2:ncol(dat)],opt.l=TRUE),
                GLM=glm(formula,family=binomial,data=dat)
  )               
  # PREDICT
  pred=as.data.frame(predict(model,newdata=alldata[,2:ncol(alldata)],type='response'),row.names=1:nrow(alldata))
  
  TSS=sensitivity(pred[,1],labels)$measure +(specificity(pred[,1],labels)$measure)-1
  TSSvalues=max(TSS)[1]
  ROCcurve=roc(pred[,1],labels)
  threshold=ROCcurve$cutoff[which(TSS==TSSvalues)[1]] ##CHECK WHAT CUTOFF VALUE IS RETURNED AND IF IT'S CORRECT
  pred=pmin(1,pred[,1]*0.5/threshold) 
  
  # CALCULATE ENTROPY OF  ALL UNLABELLED POINTS  
  E=-sum(pred*log(pred))
  return(E)
} 

# Maximize Conditional Mutual Information algorithm. It choose the point that minimize the entropy of all unlabelled instances. I.e. it minimizes prediction uncertainty based on present model
minerr_samp<-function(prior_data,predictors,formula='y~x',modeltype,labels,npoints,range) {
  library(foreach)
  # MINIMIZE ERROR THROUGH MINIMIZING ENTROPY OF UNLABELLED SET. "MCMI[MIN]" ALGORITHM FROM GUO&GREINER 2007 Optimistic Active Learning using Mutual Information
  ## THIS APPROACH SAMPLE THE POINTS WHICH WHEN ADDED TO THE TRAINING SET GIVES PREDICTION WITH THE HIGHEST ENTROPY SUMMED ACROSS ALL UNLABELED POINTS. 
  ## THIS MODEL IS THAT ONE WITH THE WIDEST SEPARATION (MORE EXTREME DIVERGENCE IN PROBABILITIES IN THE LANDSCAPE). 
  # 1) fit model adding unlabelled instances
  # 2) predict and calculate summed entropy of unlabelled instances (Sum(p*log(p))). Note that p needs to be corrected for prevalence. 
  # 3) choose unlabelled instance which minimize quantity in 2. 
  xmax=ncol(range)
  ymax=nrow(range)
  
  ind=matrix(nrow=npoints) # initialize index vector of points to sample
  data=prior_data$data # extract currently sampled data in geographic matrix format
  coords=prior_data$coords ## extract currently unsampled locations
  alldata=as.data.frame(cbind(labels[coords],predictors[coords,])) # create predictor tables and response for all unsampled data
  names(alldata)=c('PA',names(predictors))
  
 
  for (n in 1:npoints) {
  
    #E=matrix(0,length(coords)) # initialize vector of entropy values by adding each single unsampled point
    # for (i in 1:length(coords)){ # NEED A SMART WAY TO LIMIT THE NUMBER OF POINTS TO TEST SEE GUO&GREINER 2007
    #   print(i)
    #E[i]=calc_entropy_unlabelled(dataset,alldata,modeltype,threshold)
    # }
    
    library(doSNOW)
    cluster = makeCluster(8, type = "SOCK")
    registerDoSNOW(cluster)
    E <-foreach(i=1:length(coords), .combine='c',.packages=c('GRaF','AUC'),.export=c("calc_entropy_unlabelled")) %dopar% 
      calc_entropy_unlabelled(prior_data$samples[3:ncol(prior_data$samples)],formula,alldata,i,modeltype,labels) 
    stopCluster(cluster) 
    
    E=matrix(E,nrow=length(coords))
    # ORDER BY INCREASING ENTROPY. SELECT THE POINT WHICH, WHEN ADDED TO THE TRAINING DATASET, GIVES THE MODEL WHICH WHEN PROJECTED GIVES THE LOWEST SUMMED ENTROPY ACROSS ALL UNLABELED POINTS
    index=order(E)
    ind[n]=coords[index[1]]
    elim=which(coords %in% coords[index[1]])
    coords=coords[-elim]
    toadd=cbind(max(prior_data$samples$time)+1,coords[index[1]],landscape[coords[index[1]]],predictors[coords[index[1]],])
    names(toadd)=names(prior_data$samples)
    prior_data$samples=rbind(prior_data$samples,toadd)
  }
  # QUERY POINTS AND RETURN RESULTS
  data[ind]=range[ind]
  return(data)
}


### DEFINE SAMPLING FUNCTION
sampling<-function(method='uniform',time_stamp=1,prior_data,predictors,modeltype,formula,labels,landscape,
                   npoints,predictions,threshold,...){
  prevalence=sum(landscape)/length(landscape)
  # FUNCTION TAKES:
  #coords: the locations available for sampling. presently its the index of a matrix of locations
  #data: the currently available data. if none is availabe you need to specify data=matrix(NA,xmax,ymax)
  #landscape: the distribution from which you are sampling. for now this is supposed to be known
  coords=prior_data$coords
  coords_ref=coords
  data=prior_data$data
  samples=prior_data$samples
  names(samples)<-names(prior_data$samples)
  
  # FUNCTION RETURNS RESULTS LIST which is the list of new
  # 1) samples in matrix form: time in which they are sampled, the location, and the pres/abs 
  # 2) samples in geographic space
  # 3) coordinates (id) of unsampled locations
  
  sampled=matrix(nrow=0,ncol=1)
  
  ## SAMPLE POINTS DISTRIBUTION
  print('sampling')
  data<-switch(method,
               uniform=uniform_samp(xmax,ymax,coords,npoints,data,landscape),
               maxent=maxent_samp(coords,prevalence,predictions,threshold=threshold,npoints,data,landscape),
               uncert_pred=uncert_samp(coords,predictions,npoints,data,landscape),
               MIMC=minerr_samp(prior_data,predictors,formula,modeltype,labels,npoints,landscape),
               strat=stratified_samp(coords,npoints,clusters,data,landscape)
  ) 
  
  samp=as.matrix(which(data==0 | data==1))# LIST LOCATIONS THAT HAVE BEEN SAMPLED SO FAR
  print(samp)
  sampled=samp[which(!samp %in% samples$location)] # give me the new location sampled
#    print(length(samp))
#    print(length(sampled))
#    print(nrow(landscape[landscape]))
#    print(nrow(predictors[sampled,]))
  toadd=data.frame(matrix(time_stamp,npoints,1),sampled,landscape[sampled],predictors[sampled,]) # for these give me the time in which they are sampled, the location, and the pres/abs
  names(toadd)=names(samples)
  samples=rbind(samples,toadd)
  ind=which(coords_ref %in% samples$location)
  coords=coords_ref[-ind] # ELIMINATE THEM FROM THE LIST OF AVAILABLE FOR SAMPLING
  results=list(samples=samples,data=data,coords=coords)
  return(results)
}

########### FUNCTION FOR ITERATIVELY SAMPLING, FITTING AND EVALUATING ACCURACY ######################################
active_learning<-function(nsamples,method,modeltype,formula='y~x',labels,prior_data,predictors,landscape,npoints,predictions=matrix(,prod(dim(landscape)))){
  threshold=matrix(ncol=nsamples+1)
  TSSvalues=matrix(ncol=nsamples+1)
  models=vector("list",nsamples)
  AUCvalues=matrix(ncol=nsamples)
  start_time=prior_data$samples$time[nrow(prior_data$samples)]
  print('active learning')
  print(paste('sampling method =',method,'|||','modelling type =',modeltype,sep=' '))
  for (iteration in 1:nsamples){
    print(paste('iteration',iteration,sep=' '))
    pres=sum(prior_data$samples$PA)
    abs=length(prior_data$samples$PA)-pres
    print(paste('pres =',sum(prior_data$samples$PA),'abs =',abs,sep=' '))
    if (iteration==1){
      xcoord=ceiling(prior_data$samples$location/xmax)
      ycoord=(prior_data$samples$location-xmax*(xcoord-1))
      dataset=prior_data$samples
      ### FIT THE MODEL
      print('modelling for maxent sampling')
      model<-switch(modeltype,
                    GRAF=graf(dataset$PA,dataset[,4:5],opt.l=TRUE),
                    GLM=glm(formula,family=binomial,data=dataset),
                    CT=tree(as.factor(PA)~-1+temp+prec, data=dataset, split="deviance",wts=FALSE) 
      )               
      ### PREDICT OVER THE ENTIRE LANDSCAPE
      if (modeltype=='CT'){
        predictions=predict(CT,predictors,type="class")
        TSS=sensitivity(predictions,labels)$measure +(specificity(predictions,labels)$measure)-1
        TSSvalues[iteration]=max(TSS)[1]
        ROCcurve=NA#roc(predictions,labels)
      } else {
        #### COMPARE WITH "REAL" DISTRIBUTION
        predictions=as.data.frame(predict(model,newdata=predictors,type='response'),row.names=1:prod(dim(landscape)))
        TSS=sensitivity(predictions[,1],labels)$measure +(specificity(predictions[,1],labels)$measure)-1
        TSSvalues[1]=max(TSS)[1]
        ROCcurve=roc(predictions[,1],labels)
        AUCvalues[1]=auc(ROCcurve)
        threshold[1]=ROCcurve$cutoff[which(TSS==TSSvalues[iteration])[1]] ##CHECK WHAT CUTOFF VALUE IS RETURNED AND IF IT'S CORRECT
      } # END IF MODEL=CT
    } else {# IF IT>1
      ### SAMPLE  (update prior_data)
      prior_data=sampling(method,time_stamp=iteration+start_time,prior_data,predictors,modeltype,formula,labels,landscape,npoints,predictions,threshold[max(1,iteration-1)])
      
      ## RECONSTRUCT DATASET FOR MODEL FITTING FROM SAMPLES
      xcoord=ceiling(prior_data$samples$location/xmax)
      ycoord=(prior_data$samples$location-xmax*(xcoord-1))
      dataset=cbind(prior_data$samples,xcoord,ycoord)
      
      ### FIT THE MODEL
      print ('models after sampling')
      model<-switch(modeltype,
                    GRAF=graf(dataset$PA,dataset[,4:5],opt.l=TRUE),
                    GLM=glm(formula,family=binomial,data=dataset),
                    #control=tree.control(nobs=44,mincut=2,minsize=4)
                    CT=tree(as.factor(PA)~-1+temp+prec, data=dataset,control=control, split="deviance",wts=FALSE) 
                    #CT2=rpart(PA~-1+temp+prec, data=dataset, method="class", parms=list(split="gini"), control=list(minlist=3, minbucket=2))
      )               
      
      ### PREDICT OVER THE ENTIRE LANDSCAPE
      
      if (modeltype=='CT'){
        predictions=predict(CT,predictors,type="class")
        TSS=sensitivity(predictions,labels)$measure +(specificity(predictions,labels)$measure)-1
        TSSvalues[iteration]=max(TSS)[1]
        ROCcurve=NA
      } else {
        predictions=as.data.frame(predict(model,newdata=predictors,type='response'),row.names=1:prod(dim(landscape)))
        ROCcurve=roc(predictions[,1],labels)
        AUCvalues[iteration]=auc(ROCcurve)
        TSS=sensitivity(predictions[,1],labels)$measure +(specificity(predictions[,1],labels)$measure)-1
        TSSvalues[iteration]=max(TSS)[1]
        threshold[iteration]=ROCcurve$cutoff[which(TSS==TSSvalues[iteration])[1]] ##CHECK WHAT CUTOFF VALUE IS RETURNED AND IF IT'S CORRECT
        print(paste('threshold = ',threshold[iteration]))
      }
    } # end if T==1
    models[[iteration]]=model ## STORE ALL INFO ON THE BEST MODEL HER
  }# end for loop
  results=list(prior_data=prior_data,models=models,AUC=AUCvalues,TSS=TSSvalues,threshold=threshold) ## CHECK IF THE THRESHOLD WORKS
  return(results)
} # end function

