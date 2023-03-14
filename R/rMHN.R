#' rMHN(n=10000,alpha=32,beta=20,gamma=-3)
#'
#' @examples
#' ss=rMHN(n=10000,alpha=32,beta=20,gamma=-3)
#' ss=rMHN(n=1000, alpha=3, beta=1, gamma=1)
#' ss=rMHN(n=1000, alpha=3, beta=1, gamma=-1)
#' hist(ss$sample,prob=TRUE, col="grey")
#' lines(density(ss$sample))
#' @export
rMHN<-function(n=1,alpha,beta,gamma){

  # browser()
  if(gamma<=0){
    #sample=replicate(N,rMHN_negative_b(alpha=alpha,beta=beta,gamma=(gamma)))
    sample=rMHN_negative_b(n = n,  alpha=alpha,beta=beta,gamma=gamma)
    # AC_Rate=Acceptence_Rate_Neg(alpha,beta,gamma)
  }
  else if(gamma>0){
    if(alpha>1){

      #Acceptence_Rate_G=Acceptence_Rate_G(alpha,beta,gamma)
      log_K1_val_N=log_K1(alpha=alpha, beta = beta, gamma=gamma)
      log_K2_val_G=log_K2(alpha=alpha, beta = beta, gamma=gamma)
      #Acceptence_Rate_N=Acceptence_Rate_N(alpha,beta,gamma)


      if(log_K2_val_G<=log_K1_val_N){
        sample= rMHN_gamma_positive_G(n = n,alpha = alpha,beta = beta,gamma=gamma)
      }
      else if(log_K2_val_G>log_K1_val_N) {
        sample=rMHN_gamma_positive_N(n=n ,alpha = alpha,beta=beta,gamma=gamma)
        #AC_Rate=Acceptence_Rate_N
      }

    }


    if( alpha<=1){
      sample=replicate(n,rExtendedGamma(alpha,a=beta,b=(gamma)))

    }
  }
  x=sample
  return(x)
}




#' Gamma_Positive_G
#'
#' @examples
#' ss=Gamma_Positive_G(N=1000, alpha=20, beta=1, gamma=1)
#' hist(ss$sample,prob=TRUE, col="grey")
#' lines(density(ss$sample))
#' @export
rMHN_gamma_positive_G<-function(n,alpha,beta,gamma){ # N: sample size
  N=n
  # browser()
  #delta_h=beta+(2*gamma^2-sqrt(4*gamma^4+32*alpha*beta*gamma^2))/(8*alpha)
  delta_h=beta+(gamma^2-sqrt(gamma^4+8*alpha*beta*gamma^2))/(4*alpha)
  #browser()
  t <- rgamma(N,shape=alpha/2,rate=delta_h)
  x=sqrt(t)
  log_U <- log(runif(N,0,1))
  log_Rejection_Rate <- (-(beta-delta_h)*x^2+gamma*x-gamma^2/(4*(beta-delta_h)))
  sample_x <- x[log_U<=log_Rejection_Rate]  # in this step we can't gurantee the sample size =N


  size=N-length(sample_x)
  i=1
  while(size!=0){
    t <- rgamma(size,shape=alpha/2,rate=delta_h)
    x=sqrt(t)
    log_U <- log(runif(size,0,1))
    log_Rejection_Rate <-  (-(beta-delta_h)*x^2+gamma*x-gamma^2/(4*(beta-delta_h)))
    sample_x=c(sample_x,x[log_U<=log_Rejection_Rate])

    size=N-length(sample_x)

    i=i+1
    #print(i)
  }

  sample=as.vector(sample_x)
  return(sample)
}




#' Gamma_Positive_N
#'
#' @examples
#' ss=Gamma_Positive_N(N=1000, alpha=20, beta=1, gamma=1)
#' hist(ss$sample,prob=TRUE, col="grey")
#' lines(density(ss$sample))
#' @export
rMHN_gamma_positive_N<-function(n,alpha,beta,gamma){ # N: sample size
  N=n
  m=(gamma+sqrt(gamma^2+8*beta*(alpha-1)))/(4*beta);
  #m # mode

  x <- Gen_Positive_normal(n=N,mean=m, sd=1/sqrt(2*beta))
  log_U <- log(runif(N,0,1))


  log_Rejection_Rate <- (alpha-1)*log(x/m)   +((2*beta*m-gamma)*(m-x))
  condition=(log_U<=log_Rejection_Rate)
  sample <- x[condition]  # in this step we can't gurantee the sample size =N
  #browser()

  size=N-length(sample)
  i=1
  while(size!=0){
    x <- Gen_Positive_normal(n=size,mean=m, sd=1/sqrt(2*beta))
    log_U <- log(runif(size,0,1))
    log_Rejection_Rate <- (alpha-1)*log(x/m)   +((2*beta*m-gamma)*(m-x))
    condition=(log_U<=log_Rejection_Rate)
    sample0 <- x[condition]
    sample=c(sample,sample0)
    size=N-length(sample)

    i=i+1
    #print(i)
  }
  sample=as.vector(sample)
  return(sample)
}



Gen_Positive_normal<-function(n, mean, sd){
  x=rnorm(n=n, mean=mean, sd=sd)
  sample <- x[x>0]  # in this step we can't gurantee the sample size
  size=n-length(sample)

  while(size!=0){
    x=rnorm(n=size, mean=mean, sd=sd)
    sample <-c( sample,  x[x>0])
    size=n-length(sample)
    #print("in")
  }
  return(sample)
}


log_K1<-function(alpha, beta, gamma){
  mu= ( gamma + sqrt(gamma^2+8*(alpha-1)*beta)  )/(4*beta)

  log_K1_val= log(2*sqrt(pi)) +
    (alpha-1)*(  log( sqrt(beta)*(alpha-1) )  - log(2*beta*mu -gamma)   ) -(alpha-1)+ beta*mu^2

  return(log_K1_val)
}

log_K2<-function(alpha, beta, gamma){
  delta=beta+  (gamma^2-gamma*sqrt(gamma^2+8*alpha*beta)   )/(  4*alpha )

  log_K2_val= (alpha/2)*log(beta)+lgamma(alpha/2) +gamma^2/(4*(beta-delta))- (alpha/2)*log(delta)

  return(log_K2_val)
}



##################### the following function is less efficient ##################
Gamma_Positive_N_1<-function(N, alpha, beta, gamma){
  if(gamma<=0){ print("gamma must be positive to use this function."); return(NULL)}
  m=(gamma+sqrt(gamma^2+8*beta*(alpha-1)))/(4*beta);

  Gamma_Positive_N_single<-function(alp1, bet1, gam1){

    Flag=1
    while(Flag==1){
      x <- rnorm(1,mean=m, sd=1/sqrt(2*bet1))
      while(x<=0){
        x <- rnorm(1,mean=m, sd=1/sqrt(2*bet1))
      }
      log_U <- log(runif(1,0,1))
      log_Rejection_Rate <- (alp1-1)*log(x/m)   +((2*bet1*m-gam1)*(m-x))
      if(log_U<=log_Rejection_Rate){Flag=0; sample <- x}
    }
    return(sample)
  }
  val=replicate(n=N, Gamma_Positive_N_single(alp1=alpha, bet1=beta, gam1=gamma))
  return(val)
}

##################################################################################
###################################################################################


Acceptence_Rate_G<-function(alpha,beta,gamma){

  delta_h=beta+(2*gamma^2-sqrt(4*gamma^4+32*alpha*beta*gamma^2))/(8*alpha)

  set.seed(100)
  t=runif(50000,0,1000)

  fx=t^(alpha/2-1)*exp(-beta*t+gamma*sqrt(t))
  Int_fx=mean(fx)

  gx=t^(alpha/2-1)*exp(-delta_h*t)
  Int_gx=mean(gx)


  Acceptance_Rate=Int_fx/(exp(gamma^2/(4*(beta-delta_h)))*Int_gx)

  return(Acceptance_Rate)
}




Acceptence_Rate_N<-function(alpha,beta,gamma){
  m=(gamma+sqrt(gamma^2+8*beta*(alpha-1)))/(4*beta);m # mode
  set.seed(100)
  x=runif(10000,0,100)
  gx=exp(-beta*(x-m)^2)
  Int_gx=mean(gx)

  fx=x^(alpha-1)*exp(-beta*x^2+gamma*x)
  Int_fx=mean(fx)

  Acceptance_Rate=Int_fx/(m^(alpha-1)*exp(-beta*m^2+gamma*m)*Int_gx)

  return(Acceptance_Rate)
}





#ss=rMHN(N=10000,alpha=32,beta=20,gamma=-3)
#ss$sample
#ss$Acceptence_Rate
