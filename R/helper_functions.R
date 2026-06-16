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
  # Inverse gamma case
  if(a == 0){
    if(p >= 0){
      stop("Invalid parameters for GIG. When a=0, p must be negative")
    }
    mu <- -b/(2*(p+1))
  }
  # Gamma case
  if(b == 0){
    if(p <= 0){
      stop("Invalid parameters for GIG. When b=0, p must be positive")
    }
    mu <- 2*p/a
  }
  if(a != 0 & b != 0){
    mu <- sqrt(b/a)*besselK(sqrt(a*b), p+1)/besselK(sqrt(a*b), p)
  }
  return(mu)
}
rgig2.vec_fh <- Vectorize(rgig2, c("p", "a"))
var_gig <- function(p, a, b){
  # Inverse gamma case
  if(a == 0){
    if(p >= 0){
      stop("Invalid parameters for GIG. When a=0, p must be negative")
    }
    sig2 <- -b^2/(4*(p+1)^2*(p+2))
  }
  # Gamma case
  if(b == 0){
    if(p <= 0){
      stop("Invalid parameters for GIG. When b=0, p must be positive")
    }
    sig2 <- 4*p/a^2
  }
  if(a != 0 & b != 0){
    sig2 = (b/a)*(besselK(sqrt(a*b), p+2)/besselK(sqrt(a*b), p) - (besselK(sqrt(a*b), p+1)/besselK(sqrt(a*b), p))^2)
  }
  return(sig2)
}



#' @rdname gbass2bass
#' @export
gm2bm <- gbass2bass

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
rgbp<-function(v_curr,a, b, p, prop_sigma,lb, w_curr,c,r,beta_curr){
  v_cand <- exp(log(v_curr) + rnorm(1, 0, prop_sigma))
  if(v_cand >= lb){
    alpha_v <- (p*a - 1/2)*(log(v_cand) - log(v_curr)) -
      r^2/(2*w_curr*c)*(1/v_cand - 1/v_curr) -
      (a + b)*(log(1+v_cand^p) - log(1+v_curr^p)) -
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

## makes negative values 0
pos<-function(vec){
  #replace(vec,vec<0,0)
  (abs(vec)+vec)/2
}

## largest value of basis function, assuming x's in [0,1], used for scaling
const<-function(signs,knots,degree){
  cc<-prod(((signs+1)/2 - signs*knots))^degree
  if(cc==0)
    return(1)
  return(cc)
} # since a product, can find for functional & categorical pieces separately, take product

## make basis function (from continuous variables)
makeBasis<-function(signs,vars,knots,datat,degree){
  cc<-const(signs,knots,degree)
  temp1<-pos(signs*(datat[vars,,drop=F]-knots))^degree # this only works for t(data)...
  if(length(vars)==1){
    return(c(temp1)/cc)
  } else{
    temp2<-1
    for(pp in 1:length(vars)){ # faster than apply
      temp2<-temp2*temp1[pp,]
    }
    return(temp2/cc)
  }
}



