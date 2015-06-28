


plot.new()
#text(0,0.3,"Model accuracy for different sampling strategy:",cex=1.8,pos=4,xpd=T)

par(cex.lab=1.5)
par(cex.axis=1.3)
par(mgp=c(4,1,0))
#par(usr=c(0,10,0.3,0.85))
par(mar=c(5.5,6.5,3,3))

x.vals<-1:20

index=which(TSSma[1,]>0.70 & TSSma[1,]<0.90)
index=1:10
TSSma2=TSSma_graf[,index]
TSSun2=TSSunif_graf[,index]

TSSun2.l<-apply(TSSun2,1,FUN=function(x){return(quantile(x,(1/6)))})
TSSun2.m<-apply(TSSun2,1,FUN=function(x){return(quantile(x,0.5))})
TSSun2.u<-apply(TSSun2,1,FUN=function(x){return(quantile(x,(4/6)))})

TSSma2.l<-apply(TSSma2,1,FUN=function(x){return(quantile(x,(1/6)))})
TSSma2.m<-apply(TSSma2,1,FUN=function(x){return(quantile(x,0.5))})
TSSma2.u<-apply(TSSma2,1,FUN=function(x){return(quantile(x,(4/6)))})

par(las=1)
plot(x.vals,TSSma2.m,type="l",
     col="dark red",lwd=3,lty=1,ylim=c(0.3, 1),xlab=NA,ylab=NA,main='Polynomial(2), prevalence=0.2')
lines(x.vals,TSSun2.m,xlab=NA, ylab=NA, type="l",col="dark blue",lwd=3,lty=1,xaxt="n")  
par(las=0)
mtext(side = 1, text = "Additional samples", line = 3, cex=1.5)
mtext(side = 2, text = "TSS", line = 4,cex=1.5)

#axis(2,at=0.65:0:95)
#axis(1,labels=FALSE)


TSSun2.l<-apply(TSSun2,1,FUN=function(x){return(quantile(x,(1/6)))})
TSSun2.u<-apply(TSSun2,1,FUN=function(x){return(quantile(x,(4/6)))})

TSSma2.l<-apply(TSSma2,1,FUN=function(x){return(quantile(x,(1/6)))})
TSSma2.u<-apply(TSSma2,1,FUN=function(x){return(quantile(x,(4/6)))})

lines(x.vals,TSSun2.l,type="l",lty=2,col="dark blue",lwd=3)
lines(x.vals,TSSun2.u,type="l",lty=2,col="dark blue",lwd=3)

lines(x.vals,TSSma2.l,type="l",lty=2,col="dark red",lwd=3)
lines(x.vals,TSSma2.u,type="l",lty=2,col="dark red",lwd=3)

segments(x0=1,x1=1.3,y0=0.95,y1=0.95,col='dark red',lty=1,lwd=3)
segments(x0=1,x1=1.3,y0=0.91,y1=0.91,col='dark blue',lty=1,lwd=3)
text(1.3,0.95,"Uncertainty sampling",cex=1,pos=4,xpd=T)
text(1.3,0.91,"Uniform random sampling",cex=1,pos=4,xpd=T)















a=predict(results_prior$models[[1]],newdata=predictors,type='response')
a=signif(matrix(a[,1],nrow=15,ncol=15),digits=4)
rbPal <- colorRampPalette(c('light blue','blue', 'red','dark red'))
col <- rbPal(45)[as.numeric(cut(a,breaks = 45))]
#plot(covar[,1],-covar[,2],type='p',pch=22, col=col,bg=col,cex=4)

plot(predictors[,2],-predictors[,1],xlab='Predictor 1',ylab='Predictor 2',type='p',pch=22, col=col,bg=col,cex=4)


par(mfrow=c(1,3))
plot3d(results_prior$models[[1]], dims = c(1, 2), prior = FALSE, peak = TRUE)

### PRINT PREDICTIONS THIS WILL BE IN GEOGRAPHIC SPACE WITH SP AND RGEOS
a=predict(results_maxent[[2]][[10]],newdata=predictors,type='response')
a=signif(matrix(a[,1],nrow=15,ncol=15),digits=4)
rbPal <- colorRampPalette(c('red','blue'))
col <- rbPal(15)[as.numeric(cut(a,breaks = 15))]
plot(covar[,1],-covar[,2],type='p',pch=22, col=col,bg=col,cex=4)

## PLOT PREDICTIONS (GEOGRAPHIC SPACE)
plot.new()
# axis(1)
# axis(2)
par(mfrow=c(1,1))
plot(predictors[,1],-predictors[,2],type='p',pch=22, col='black',bg=col,cex=4)

# par(mar=c(5, 5, 5, 5))
# par(pin=c(10, 10))
# plot.window(xlim=c(0,15), ylim=c(-15,0))



## PLOT RESPONSES (ENVIRONMENTAL SPACE)
par(mfrow = c(1, 2))
plot(results_uniform[[2]][[4]])
plot(results_maxent[[2]][[4]])
plot(results_prior[[2]][[1]])
meant=mean(predictors[,1])
meanp=mean(predictors[,2])
probs_temp=1/(1+exp(-(interc+t1*(1:15)+t2*(1:15)^2+t3*(1:15)^3+
                   p1*meanp+p2*meanp^2+p3*meanp^3)))
probs_prep=1/(1+exp(-(interc+t1*meant+t2*meant^2+t3*meant^3+
                        p1*(1:15)+p2*(1:15)^2+p3*(1:15)^3)))
lines(1:15,probs_temp)
lines(1:15,probs_prep)

par(mfrow = c(2, 6))
par(cex = 0.6)
par(mar = c(3, 3, 0, 0), oma = c(1, 1, 1, 1))
for (i in c(1,5,9,13,17,20)) {
  if (i==1){
    plot(results_prior[[2]][[1]])
    mtext(i, side = 3, line = -1, adj = 0.1, cex = 0.6)
  }else {
    plot(results_uniform[[2]][[i]])
    mtext(i, side = 3, line = -1, adj = 0.1, cex = 0.6)
  }
}


par(mfrow = c(2, 6))
par(cex = 0.6)
par(mar = c(3, 3, 0, 0), oma = c(1, 1, 1, 1))
for (i in c(1,5,9,13,17,20)) {
  if (i==1){
    plot(results_prior[[2]][[1]])
    mtext(i, side = 3, line = -1, adj = 0.1, cex = 0.6)
  }else {
    plot(results_maxent[[2]][[i]])
    mtext(i, side = 3, line = -1, adj = 0.1, cex = 0.6)
  }
}

real_prec=1/(1+exp(-(interc+t1*mean(predictors[,1])+t2*mean(predictors[,1])^2+t3*mean(predictors[,1])^3+
                       p1*(1:15)+p2*(1:15)^2+p3*(1:15))))