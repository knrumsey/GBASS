#' Generalized Inverse Gaussian Generator
#'
#' This function generates samples from the GIG(p, a, b) distribution with density function
#' f(x) = cx^(p-1)exp(-1/2*(a*x + b/x))
#'
#' @param p a real valued parameter
#' @param a a non-negative parameter. If a=0, then p must be negative.
#' @param b a non-negative parameter. If b=0, then p must be positive.
#' @details A uniformly bounded rejection sample based on Hörmann et. al. (2014). Special cases include Gamma (b=0, p>0), Inverse Gamma (a=0, p<0) and Inverse Gaussian (p=-1/2).
#' @references Hörmann, Wolfgang, and Josef Leydold. "Generating generalized inverse Gaussian random variates." Statistics and Computing 24.4 (2014): 547-557.
#' @export
#' @examples
#' x <- rep(NA, 1000)
#' for(i in 1:1000) x[i] <- rgig2(2, 3, 0.5)
#' hist(x)
#'
rgig2 <- function(p, a, b){
  res <- ifelse(a==0, 1/rgamma(1, -p, b/2), GIGrvg::rgig(1, lambda=p, psi=a, chi=b))
}
rgig2.vec <- Vectorize(rgig2, "b")
mu_gig <- function(p, a, b){
  sqrt(b/a)*besselK(sqrt(a*b), p+1)/besselK(sqrt(a*b), p)
}
var_gig <- function(p, a, b){
  (b/a)*(besselK(sqrt(a*b), p+2)/besselK(sqrt(a*b), p) - (besselK(sqrt(a*b), p+1)/besselK(sqrt(a*b), p))^2)
}

#' Convert GBASS to BASS
#'
#' This function converts an object of class `gbass` to an object of class `bass`. This is convenient for taking advantage of
#' the sobol decomposition functions available in the BASS package.
#'
#' @param gm an object of class gbass
#' @export
#' @examples
#' #The following are equivalent
#' n <- 100 #Number of observations
#' p <- 4   #Number of variables (beyond p = 2, variables are inert)
#' X <- matrix(runif(n*p), nrow=n)
#' y <- apply(X, 1, ff1)
#' gm <- gbass(X, Y, nmcmc=1000, nburn=901)
#' bm <- gm2bm(gm)
#' sob <- sobol(bm)
#' plot(sob)
#'
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

  out$xx.des <- gm$X
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


dmwnchBass<-function(z.vec,vars){
  z.vec <- z.vec/sum(z.vec)
  alpha<-z.vec[vars]/sum(z.vec[-vars])
  j<-length(alpha)
  ss<-1 + (-1)^j * 1/(sum(alpha)+1)
  if(j > 1){
    for(i in 1:(j-1))
      ss <- ss + (-1)^(i) * sum(1/(colSums(combn(alpha,i))+1))
  }
  ss
}

logdet <- function(A) Matrix::determinant(A)$modulus
move_type <- function(M, Mmax, Pmove){
  if(M <= 0)    return("B")
  if(M >= Mmax) return(sample(c("D", "M"), 1, prob=Pmove[2:3]))
  return(sample(c("B", "D", "M"), 1, prob=Pmove))
}
rgbp<-function(v_curr,v_prior,w_curr,c,r,beta_curr){
  v_cand <- exp(log(v_curr) + rnorm(1, 0, v_prior$prop_sigma))
  if(v_cand >= v_prior$lb){
    alpha_v <- (v_prior$p*v_prior$a - 1/2)*(log(v_cand) - log(v_curr)) -
      r^2/(2*w_curr*c)*(1/v_cand - 1/v_curr) -
      (v_prior$a + v_prior$b)*(log(1+v_cand^v_prior$p) - log(1+v_curr^v_prior$p)) -
      beta_curr^2/(2*c)*(v_cand - v_curr)
  }else{
    alpha_v <- -Inf
  }
  if(log(runif(1)) < alpha_v){
    return(v_cand)
  }
  return(v_curr)
}
rgbp.vec<-Vectorize(rgbp,c("v_curr","r"))
TCP  <- function(...) Matrix::tcrossprod(...)
symchol <- function(...) base::chol(...)

myTimestamp <-function(){
  x<-Sys.time()
  paste('#--',format(x,"%b %d %X"),'--#')
}




