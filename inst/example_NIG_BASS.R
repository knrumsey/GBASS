library(lhs)
tform <- function(x){
  a <- c(.05, 100, 63070,990,63.1,700,1120,9855)
  b <- c(.15, 50000,115600, 1110, 116, 820, 1680, 12045)
  x*(b-a) + a
}

borehole <- function(xx)
{
  ##########################################################################
  #
  # BOREHOLE FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it
  # and/or modify it under the terms of the GNU General Public License as
  # published by the Free Software Foundation; version 2.0 of the License.
  # Accordingly, this program is distributed in the hope that it will be
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # OUTPUT AND INPUT:
  #
  # y  = water flow rate
  # xx = c(rw, r, Tu, Hu, Tl, Hl, L, Kw)
  #
  ##########################################################################

  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]

  frac1 <- 2 * pi * Tu * (Hu-Hl)

  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)

  y <- frac1 / frac2
  return(y)
}

# rw ∈ [0.05, 0.15] 	radius of borehole (m)
# r ∈ [100, 50000] 	radius of influence (m)
# Tu ∈ [63070, 115600]    	transmissivity of upper aquifer (m2/yr)
# Hu ∈ [990, 1110] 	potentiometric head of upper aquifer (m)
# Tl ∈ [63.1, 116] 	transmissivity of lower aquifer (m2/yr)
# Hl ∈ [700, 820] 	potentiometric head of lower aquifer (m)
# L ∈ [1120, 1680] 	length of borehole (m)
# Kw ∈ [9855, 12045] 	hydraulic conductivity of borehole (m/yr)
rr<-matrix(c(
  0.05, 0.15,
  100, 50000,
  63070, 115600,
  990, 1110,
  63.1, 116,
  700, 820,
  1120, 1680,
  9855, 12045
),nrow=2)

library(lhs)
x<-lhs::randomLHS(200,8)
x.un<-do.call(cbind,lapply(1:8,function(i) BASS:::unscale.range(x[,i],rr[,i])))

y<-apply(x.un,1,borehole)

xx<-randomLHS(1000,8)
xx.un<-do.call(cbind,lapply(1:8,function(i) BASS:::unscale.range(xx[,i],rr[,i])))
yy<-apply(xx.un,1,borehole)

library(BASS)
bmod<-bass(x,y,h1=1,h2=.001)
plot(bmod)
library(GBASS)
library(Matrix)
gmod3<-gbass(x,y,v_prior = build_prior('GIG',c(5,10,0),.5),w_prior = build_prior('GIG',c(0,0,0),.5))
gmod<-gbass(x,y,v_prior = build_prior('GIG',c(-15,0,30),.5),w_prior = build_prior('GIG',c(0,0,0),.5)) # vi ~ IG(15,15)
sqrt(mean((yy-rowMeans(predict(gmod,xx)))^2))
sqrt(mean((yy-rowMeans(predict(gmod3,xx)))^2))
sqrt(mean((yy-colMeans(predict(bmod,xx)))^2))

library(parallel)
pp<-build_prior('GIG',c(0,0,0),.5)
pp$prop_sigma=.5
gm.list<-mclapply(c(.1,.25,.5,.75,.9),function(q) qbass(x,y,q,w_prior = pp),mc.cores=3)


gm2bm<-function(gm){
  out<-list()
  out$degree<-1

  nmcmc<-length(gm$a)
  maxb<-max(gm$M)
  beta<-matrix(nrow=nmcmc,ncol=maxb+1)
  for(i in 1:nmcmc)
    beta[i,1:length(gm$a[[i]])]<-gm$a[[i]]

  out$beta<-beta
  out$nbasis<-gm$M
  out$p<-8
  gm$lookup[[1]]

  model.lookup<-1
  for(i in 2:nmcmc){
    if(gm$M[i]==gm$M[i-1]){
      if(all(gm$basis[[i]] == gm$basis[[i-1]])){
        model.lookup[i]<-model.lookup[i-1]
      } else{
        model.lookup[i]<-model.lookup[i-1]+1
      }
    } else{
      model.lookup[i]<-model.lookup[i-1]+1
    }
  }
  out$model.lookup<-model.lookup
  out$n.models<-max(model.lookup)
  out$des<-T
  out$func<-F
  out$cat<-F
  max.int<-3
  n.int.des<-matrix(nrow=out$n.models,ncol=maxb)
  signs.des<-vars.des<-knots.des<-array(dim = c(out$n.models,maxb,max.int))
  for(i in 1:out$n.models){
    ind<-which(model.lookup==i)
    for(j in 1:gm$M[ind[1]]){
      n.int.des[i,j]<-gm$lookup[[gm$basis[ind[1]][[1]][j]]]$J
      signs.des[i,j,1:n.int.des[i,j]]<-gm$lookup[[gm$basis[ind[1]][[1]][j]]]$s
      vars.des[i,j,1:n.int.des[i,j]]<-gm$lookup[[gm$basis[ind[1]][[1]][j]]]$u
      knots.des[i,j,1:n.int.des[i,j]]<-gm$lookup[[gm$basis[ind[1]][[1]][j]]]$t
    }
  }

  out$n.int.des<-n.int.des
  out$signs.des<-signs.des
  out$knots.des<-knots.des
  out$vars.des<-vars.des
  out$cx<-rep('numeric',out$p)
  out$range.des<-rbind(rep(0,out$p),rep(1,out$p))
  out$nburn<-0
  out$thin<-1
  out$nmcmc<-nmcmc
  out$pdes=out$p
  out$pcat<-0
  out$pfunc<-0
  out$maxInt.des<-max.int
  out$maxInt.cat<-0
  out$maxInt.func<-0

  class(out)<-'bass'

  return(out)
}


set.seed(105)
n <- 500
X <- randomLHS(n, 8)
y <- apply(X, 1, function(z) borehole(tform(z)))
y2 <- y + 4.5*(rgamma(n, .1, .1) - 1) + rnorm(n, 0, .5)
mod1 <- bass(X, y2)
mod3 <- nwbass(X, y2, m_beta=0, s_beta=10, m_gamma=0,
               s_gamma=3, a_lambda=.1, b_lambda=.1,
               nmcmc=50000, nburn=46001, thin=4,
               w_prior=list(type="GIG", p=-2, a=0, b=2, prop_sigma=0.2))

#Test set
n <- 5000
XX <- randomLHS(n, 8)
yy <- apply(XX, 1, function(z) borehole(tform(z)))
yy2 <- yy + 4.5*(rgamma(n, .1, .1) - 1) + rnorm(n, 0, 0.5)

#CHECK COVERAGE OF BASS
cnt <- 0
wid <- rep(NA, n)
for(i in 1:n){
  yhat = predict(mod1, XX[i,])
  CI <- quantile(yhat + rnorm(1000, 0, sqrt(mod1$s2)),
                 probs=c(.05, .95))
  if(yy2[i] > CI[1] & yy2[i] < CI[2]){
    cnt <- cnt + 1/n
  }
  wid[i] <- diff(CI)
}
cat("BASS obtained coverage:", cnt, "with a median CI width of", median(wid))

#CHECK COVERAGE OF NW-BASS
library("TeachingDemos")
cnt3 <- 0
wid3 <- rep(NA, n)
yhatM = predict(mod3, (as.matrix(XX)))
for(i in 1:n){
  yhat = yhatM[i,]
  eps <- rep(NA, 1000)
  for(j in 1:1000){
    vv <- rgig2(-1/2,mod3$gam[j]^2,1)
    eps[j] <- sqrt(mod3$w[j])*(mod3$bet[j]*vv + sqrt(vv)*rnorm(1)) #- mod3$bet[j]/mod3$gam[j]
  }
  #CI <- quantile(yhat + eps, probs=c(.05, .95))
  CI <- emp.hpd(yhat+eps, conf=0.95)
  if(yy2[i] > CI[1]  & yy2[i] < CI[2]){
    cnt3 <- cnt3 + 1/n
  }
  wid3[i] <- diff(CI)
  if(i%%100 == 0) print("hey")
}
cat("NW-BASS obtained coverage:", cnt3, "with a median CI width of", median(wid3))
hist(wid/wid3, main="",
     xlab="Ratio of CI Widths", cex.lab=1.3, cex.main=1.4)



#SIMULATE DISTRIBUTION
library(RColorBrewer)
bob <- brewer.pal(3, "Set1")
eps <- rep(NA, 5000)
for(i in 1:1000){
  ii <- 5*(i-1)
  for(j in 1:5){
    vv <- rgig2(-1/2,mod3$gam[i]^2,1)
    eps[ii+j] <- sqrt(mod3$w[i])*(mod3$bet[i]*vv + sqrt(vv)*rnorm(1) - mod3$bet[i]/mod3$gam[i])
  }
}
nn <- 10000
eps_true <- 4.5*(rgamma(nn, .1, .1) - 1) + rnorm(nn, 0, 0.5)
lines(density(eps_true-mean(eps_true), bw=2, n=1001), lwd=2, col=bob[2])
eps_bass <- rnorm(10000, 0, sqrt(mod1$s2))

plot(NULL, xlim=c(-10, 30), main="", ylab="", col="white", lwd=2, ylim=c(0, .18))
grid()
polygon(density(eps_true-0*mean(eps_true),bw=2,n=5001), col=adjustcolor(bob[1], alpha.f=0.6), border='white')
polygon(density(eps-0*mean(eps),bw=2,n=5001), col=adjustcolor(bob[2], alpha.f=0.6), border='white')
polygon(density(eps_bass-0*mean(eps_bass),bw=2,n=5001), col=adjustcolor(bob[3], alpha.f=0.6), border='white')
legend('topright', c("True Residual Distribution", "Estimated Distribution (Gauss)", "Estimated Distribution (NW)"),
       lwd=3, fill=adjustcolor(bob[c(1,3,2)], alpha.f=0.8), lty=-1,
       cex=1.5, bty="n")

