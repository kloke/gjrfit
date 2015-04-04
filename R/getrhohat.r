getrhohat <-
function(x,y,block,visit,eps=0.01) { 

nvec<-tapply(!is.na(y),block,sum)
ublock<-unique(block)
m<-length(ublock)

myorder<-order(block,visit)

y<-y[myorder]
x<-x[myorder,]

#ehat<-rfitqr(x,y)$residuals
ehat<-rfit(y~x,TAU='N')$residuals

Ehat<-matrix(nrow=m,ncol=max(nvec))
for( k in 1:m ) {
	ek<-ehat[block==ublock[k]]
	vk<-visit[block==ublock[k]]
	Ehat[k,vk]<-ek
}

y<-as.vector(t(Ehat[,2:ncol(Ehat)]))
x<-as.vector(t(Ehat[,1:(ncol(Ehat)-1)]))

ind<-(!is.na(y))*(!is.na(x))
y<-y[ind==1]
x<-x[ind==1]

yhat1<-hbrfit(y~x)$fitted.values
rhohat<-lsfit(x,yhat1,intercept=FALSE)$coef

rhohat <- ifelse(abs(rhohat) >= 1, sign(rhohat)*(1-eps),rhohat)

rhohat

}
