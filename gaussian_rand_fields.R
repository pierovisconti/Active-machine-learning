dataset2=cbind(dataset[,4:5],dataset[,4]*dataset[,5])
names(dataset2)=c('xcoord','ycoord','interaction')
covar=cbind(xcoord,ycoord)
library(GRaF)  ## fit Gaussian Random Fields to spatially autocorrelated presence data
gauss_model_opt=graf(dataset$PA,dataset[,4:5],opt.l=TRUE)
gauss_model=graf(dataset$PA,dataset[,4:5],opt.l=FALSE)
print(gauss_model)
par(mfrow = c(1, 2))
plot(gauss_model)
par(mfrow = c(1, 2))
plot(gauss_model_opt)
pred_opt <- predict(gauss_model_opt, covar)
pred <- predict(gauss_model, covar)
landscape_vect=matrix(c(t(landscape)),nrow=nrow(landscape)*ncol(landscape))
cols=c('white','black')
#plot.new()
plot(covar[,1],-covar[,2],type='p',pch=22, col='black',bg=cols[landscape_vect+1],cex=4)

## COMPARE DIFFERENT REGRESSION MODELS AND THEIR PERFORMANCE
prior_data=sampling('maxent',n=1,prior_data,predictors,'GRAF',labels,landscape,4,predictions,threshold=0.5)
gauss_model_opt=graf(prior_data$samples$PA,prior_data$samples[,4:5],opt.l=TRUE)
plot(gauss_model_opt)
model=glm(PA~-1+temp+prec, family=binomial(link='logit'),data=prior_data$samples)
#predictions=as.data.frame(predict(model,newdata=predictors,type='response'),row.names=1:prod(dim(landscape)))
predictions=as.matrix(predict(model,newdata=predictors,type='response'))
predictions_gauss=as.matrix(predict(gauss_model_opt,predictors))
#plot.new()

#plot(covar[,1],-covar[,2],type='p',pch=22, col='black',bg=cols[round(predictions)+1],cex=4)
rbPal <- colorRampPalette(c('red','blue'))
Col <- rbPal(10)[as.numeric(cut(predictions_gauss[,1],breaks = 10))]
plot(covar[,1],-covar[,2],type='p',pch=22, col=Col,bg=Col,cex=4)




model=glm(PA~-1+xcoord+ycoord, family=binomial(link='probit'),data=dataset)
predictions_logit=as.data.frame(predict(model,newdata=covar,type='response'))
predictions_probit=as.data.frame(predict(model,newdata=cbind[xcoord,ycoord],type='response'))
preds=cbind(landscape_vect,pred[,1],pred_opt[,1],predictions_logit,predictions_probit)
names(preds)=c('pa','randfi','randfi.opt','logit','probit')
preds