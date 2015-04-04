gjrfit.ar1 <-
function(x,y,block,visit) { 

### NOTE: assumes that 1 is in the column space of x ###

x<-as.matrix(x)

ublock<-unique(block)
m <- length(ublock)

### estimate the AR(1) parameter ###
rhohat<-getrhohat(x,y,block,visit)

### calculate ystar & Xstar ###
ys<-vector(mode='numeric',length=length(y))
xs<-matrix(nrow=nrow(x),ncol=ncol(x))

for( k in 1:m ) {

	yk<-y[block==ublock[k]]
	xk<-x[block==ublock[k],]
	vk<-visit[block==ublock[k]]
	Shat<-rhohat^abs(outer(vk,vk,'-'))
	eigShat<-eigen(Shat)
	ShatInv<-eigShat$vectors%*%diag(1/sqrt(eigShat$values))%*%t(eigShat$vectors)
	ys[block==ublock[k]]<-as.vector(ShatInv%*%yk)
	xs[block==ublock[k],]<-ShatInv%*%xk

}

### yhat1s = alphahat 1 + xc* betahat ###
# fit<-rfitqr(xs,ys,fitint=FALSE)
fit <- rfit(ys~xs,tau='N')
yhats1<-fit$fitted.values


Qw<-qr.Q(fit$qrx1)

# ehats<-ys-yhats1
# these are the iid residuals (to be used in the bootstrap)
ehats <- fit$residuals
r <- rank(ehats, ties.method = "first")/(length(ehats) + 1)
scorehat<-getScores(fit$scores,r)

tauhat<-gettau(ehats,fit$qrx1$rank,scores=fit$scores)

### yhats = Hx* yhats1 (project onto correct space) ###
qrxs<-qr(xs)
Qxs<-qr.Q(qrxs)
yhats<-qr.fitted(qrxs,yhats1)

### transform back to correct space ###
yhat<-vector(mode='numeric',length=length(y))

for( k in 1:m ) {

	yk<-yhats[block==ublock[k]]
	vk<-visit[block==ublock[k]]
	Shat<-rhohat^abs(outer(vk,vk,'-'))
	eigShat<-eigen(Shat)
	Shat<-eigShat$vectors%*%diag(sqrt(eigShat$values))%*%t(eigShat$vectors)
	yhat[block==ublock[k]]<-as.vector(Shat%*%yk)

}

### solve ls problem x betahat = yhat ###
betahat<-lsfit(x,yhat,intercept=FALSE)$coef

### calc variance of betahat ###

xpxi<-chol2inv(chol(crossprod(x)))

A1<-matrix(0,nrow=ncol(x),ncol=ncol(Qxs))

A2<-matrix(0,nrow=ncol(Qw),ncol=ncol(Qw))

for( k in 1:m ) {

	ind<-block==ublock[k]

	vk<-visit[ind]
	xk<-x[ind,]
	Qxsk<-Qxs[ind,]
	Qwk<-Qw[ind,]
	Shat<-rhohat^abs(outer(vk,vk,'-'))
	eigShat<-eigen(Shat)
	Shat<-eigShat$vectors%*%diag(sqrt(eigShat$values))%*%t(eigShat$vectors)

	A1<-A1+t(xk)%*%Shat%*%Qxsk

	scorehatk<-scorehat[ind]
	a2<-crossprod(Qwk,scorehatk)

	A2<-A2+tcrossprod(a2)

}

B<-crossprod(Qxs,Qw)

A1B<-A1%*%B

varbetahat<-tauhat*tauhat*xpxi%*%A1B%*%A2%*%t(A1B)%*%xpxi

list(betahat=betahat,varbetahat=varbetahat,rhohat=rhohat,fitted.values=yhat,residuals=y-yhat,residuals.iid=ehats,block=block,visit=visit,x=x,y=y,DF=m)

}
