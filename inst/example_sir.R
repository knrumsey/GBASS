library(EpiModel)
library(parallel)

##############################
## get training data - 2000 model runs with 3 parameters (but only 2 active), 30 time steps
n<-2000
p<-3
nt<-21
X<-matrix(runif(n*p),ncol=p) # inputs on 0-1 scale
y<-matrix(nrow=nrow(X),ncol=nt)
#Parameter range chosen so that P(# New Infections > 0) >= 0.05
for(i in 1:n){
  param <- param.icm(inf.prob = X[i,1]*.5 + .25, act.rate = X[i,2]*1 + .5) # rescale inputs...
  init <- init.icm(s.num = 5000, i.num = 1) # 500 suceptible, 1 infected
  control <- control.icm(type = "SI", nsims = 1, nsteps = nt) # single simulation
  mod <- icm(param, init, control)
  y[i,]<-unlist(mod$epi$i.num)
}
matplot(t(y),type='l')

##############################
## get validation data
ntest<-100
nreps<-1000 # number of replications for each parameter combination
XX<-matrix(runif(ntest*p),ncol=p)
yy<-array(dim=c(nrow(XX),nt,nreps))
for(i in 1:ntest){
  param <- param.icm(inf.prob = XX[i,1]*.5+.25, act.rate = XX[i,2]*1+.5)
  init <- init.icm(s.num = 5000, i.num = 1)
  control <- control.icm(type = "SI", nsims = nreps, nsteps = nt)
  mod <- icm(param, init, control)
  yy[i,,]<-as.matrix(mod$epi$i.num)
}

##############################
## plot a single set of 10 replications, 300 time steps
param <- param.icm(inf.prob = .2, act.rate = .25)
init <- init.icm(s.num = 5000, i.num = 1)
control <- control.icm(type = "SI", nsims = 10, nsteps = 300)
mod <- icm(param, init, control)
plot(mod)


QQ_list<-c(.05,.1,.175,.25,.5,.75,.825,.9,.95)

##############################
## fit quantile models for 30th time point
library(GBASS)
pp<-build_prior('GIG',c(0,0,0),.5)
pp$prop_sigma=.5
gm.list<-mclapply(QQ_list,function(q) qbass(X,y[,21],q,w_prior = pp),mc.cores=3)

##############################
## predict validation data, compare to estimated quantiles

# not sure how to include noise
par(mfrow=c(2,3))
for(i in 1:length(QQ_list)){
  pp<-predict(gm.list[[i]],XX) #+ t(matrix(rnorm(200*1000,sd=sqrt(gm.list[[i]]$w)),nrow=1000))
  qq<-apply(pp,1,quantile,probs=c(.025,.975))
  qtrue<-apply(yy[,21,],1,quantile,probs=QQ_list[i])
  plot(qtrue,rowMeans(pp),main=QQ_list[i])
  segments(qtrue,qq[1,],qtrue,qq[2,],col='lightgrey')
  points(qtrue,rowMeans(pp)); abline(a=0,b=1,col=2)
}

##############################
## translate to bass structure, get sobol
gbm<-lapply(gm.list,gm2bm)
sob<-lapply(gbm,GBASS:::sobol)

for(i in 1:length(QQ_list)){
  #plot(sob[[i]])
  plot()
}

##############################
## compare to BASS, quantiles of posterior predictive dist
library(BASS)
bmod<-bass(X,y[,21])
plot(bmod)

# mean prediction
pp<-predict(bmod,XX)
plot(apply(yy[,21,],1,mean),colMeans(pp)); abline(a=0,b=1,col=2)

# sobol
plot(sobol(bmod))

#Make sobol plot
library(RColorBrewer)
bob <- brewer.pal(7, "Dark2")
SS <- array(NA, dim=c(length(QQ_list), 7, 3))
for(i in 1:length(QQ_list)){
  SS[i,,] <- t(apply(sob[[i]]$S, 2, quantile, prob=c(.025, .5, .975)))
}
par(mfrow=c(1,2))
plot(NULL, xlim=range(QQ_list), ylim=1.2*range(SS), xlab='Quantile', ylab='Sensitivity', cex.lab=1.2, main="Sensitivity")
grid()
for(j in 1:7){
  polygon(c(QQ_list, rev(QQ_list)), c(SS[,j,1], rev(SS[,j,3])),
          col=adjustcolor(bob[j], alpha.f=0.4), border='white')
  lines(QQ_list, SS[,j,2], lwd=2, col=bob[j])
}
#legend("topright", c("x1", "x2", "x3", "x1*x2", "x1*x3", "x2*x3", "x1*x2*x3"),
#       col=bob, lwd=2, bty='n', horiz=FALSE)
legend(0.1, 0.7, c("x1", "x2", "x3"),
       col=bob, lwd=2, bty='n', horiz=FALSE)
legend(0.5, 0.703, c("x1*x2", "x1*x3", "x2*x3", "x1*x2*x3"),
       col=bob[4:8], lwd=2, bty='n', horiz=FALSE)


library(RColorBrewer)
bob <- brewer.pal(7, "Dark2")
TT <- array(NA, dim=c(length(QQ_list), 3, 3))
for(i in 1:length(QQ_list)){
  TT[i,,] <- t(apply(sob[[i]]$T, 2, quantile, prob=c(.025, .5, .975)))
}
plot(NULL, xlim=range(QQ_list), ylim=c(0, 0.9), xlab='Quantile', ylab='Sensitivity', cex.lab=1.2, main="Total Sensitivity")
grid()
for(j in 1:3){
  polygon(c(QQ_list, rev(QQ_list)), c(TT[,j,1], rev(TT[,j,3])),
          col=adjustcolor(bob[j], alpha.f=0.4), border='white')
  lines(QQ_list, TT[,j,2], lwd=2, col=bob[j])
}
#legend("topright", c("x1", "x2", "x3", "x1*x2", "x1*x3", "x2*x3", "x1*x2*x3"),
#       col=bob, lwd=2, bty='n', horiz=FALSE)
legend("topright", c("x1", "x2", "x3"),
       col=bob, lwd=2, bty='n', horiz=FALSE)






# quantiles of posterior predictive dist compared to true quantiles
for(i in 1:length(QQ_list)){
  pp<-predict(bmod,XX)+matrix(rnorm(1000*ntest,sd = sqrt(bmod$s2)),nrow=1000)
  qq<-apply(pp,2,quantile,probs=QQ_list[i])
  qtrue<-apply(yy[,21,],1,quantile,probs=QQ_list[i])
  plot(qtrue,qq,main=QQ_list[i])
  abline(a=0,b=1,col=2)
}




sobol_q_plot <- function(X, y, q_list = c(0.05, 0.25, 0.5, 0.75, 0.9), cores=1, alpha=0.05, colors="orange", ...){

  if(cores == 1){
    gm.list <- lapply(q_list, function(q) qbass(X, y, q), ...)
  }else{
    gm.list <- mclapply(q_list, function(q) qbass(X, y, q, ...), mc.cores=cores)
  }
  ## translate to bass structure, get sobol
  gbm<-lapply(gm.list,gm2bm)
  sob<-lapply(gbm,GBASS:::sobol)

  nS <- ncol(sob[[1]]$S)
  nT <- ncol(sob[[1]]$T)

  #Gather sobol data
  SS <- array(NA, dim=c(length(q_list), nS, 3))
  TT <- array(NA, dim=c(length(q_list), nT, 3))
  for(i in seq_along(q_list)){
    SS[i,,] <- t(apply(sob[[i]]$S, 2, quantile, prob=c(alpha/2, .5, 1-alpha/2)))
    TT[i,,] <- t(apply(sob[[i]]$T, 2, quantile, prob=c(alpha/2, .5, 1-alpha/2)))
  }

  #Make sensitivity plot
  par(mfrow=c(1,1))
  plot(NULL, xlim=range(q_list), ylim=range(SS), xlab='Quantile', ylab='Sensitivity', cex.lab=1.2, main="Sensitivity")
  grid()
  for(j in 1:nS){
    polygon(c(q_list, rev(q_list)), c(SS[,j,1], rev(SS[,j,3])),
            col=colors[j], border='white')
    lines(q_list, SS[,j,2], lwd=2, col=colors[j])
  }

  #Make total sensitivity plot
  plot(NULL, xlim=range(q_list), ylim=range(TT), xlab='Quantile', ylab='Sensitivity', cex.lab=1.2, main="Total Sensitivity")
  grid()
  for(j in 1:nT){
    polygon(c(q_list, rev(q_list)), c(TT[,j,1], rev(TT[,j,3])),
            col=colors[j], border='white')
    lines(q_list, TT[,j,2], lwd=2, col=colors[j])
  }
  par(mfrow=c(1,1))

  return(sob)

}






