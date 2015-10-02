library(parallel)
library(Rcpp)
sourceCpp("num_records.cpp")
library(skewt)

options(mc.cores=detectCores())

mysumratio=function(m){
  k=1:(m-1)
  mytmp=sum(sqrt(k)/(m-k))/(2*sqrt(m)*(atanh(sqrt(1-1/m))-sqrt(1-1/m))+0.5*sqrt(m-1))
  return(mytmp)
}
mysumratioV=Vectorize(mysumratio)

myapproxZero=function(z,nmax=100000){
  n=1:nmax
  return(2/sqrt(pi)*sum( sqrt(n) *z^n ))
}

myapproxZeroV=Vectorize(myapproxZero)


myapproxSt=function(z,nmax=100000){
  n=1:nmax
  return(2/sqrt(pi)*sum( 2*(sqrt(n)*(atanh(sqrt(1-1/n))-sqrt(1-1/n)))+0.5*sqrt(n-1) *z^n ))
}

myapproxStV=Vectorize(myapproxSt)

T_moments=function(k,nu){
  return(1/(sqrt(pi)*gamma(nu/2))*gamma((k+1)/2)*gamma((nu-k)/2)*nu^(k/2))
}

exp_num_records_approx=function(n){
  return(sqrt(4*n/pi))
}

prob_num_records=function(R,n){
  K=2*n-R+1
  return(choose(K,n)*2^(-K))
}

exp_num_records=function(n){
  return(choose(2*n,n)*(2*n+1)/2^(2*n))
}

exp_sd_records=function(n){
  return(sqrt((2-4/pi)*n))
}

igamma=function(a,x) {
  return(pgamma(x,a,lower=FALSE)*gamma(a))
}

approx_diff_driftT=function(n,alpha,c,sigma){
  return(approx_diff_driftGauss(n,alpha,c,sigma)+
           4/(sqrt(3)*pi^1.5)*c/sigma*(2*sqrt(n)*(atanh(sqrt(1-1/n))-sqrt(1-1/n))))
  
}


approx_diff_driftGauss=function(n,alpha,c,sigma){
  return(c/sigma*sqrt(2)/pi*(n*atan(sqrt(n-1))-0.5*sqrt(n-1)))
}



estimateSNRFromNumRecordsStudent=function(n,numRecords){
  #   f=function(x) numRecords-sqrt(4*n/pi)-x*sqrt(2)/pi*(n*atan(sqrt(n))-sqrt(n))-x*16/(sqrt(3)*pi^1.5)*sqrt(n)*(sqrt(1+1/n)*asinh(sqrt(n))-1)
  #   return(uniroot(f,interval=c(-3,3))$root)
  
  SNR=(numRecords-sqrt(4*n/pi))/(sqrt(2)/pi*(n*atan(sqrt(n))-sqrt(n))+8/(sqrt(3)*pi^1.5)*sqrt(n)*(atanh(sqrt(1-1/n))-sqrt(1-1/n)))
  return(SNR)
}  

estimateSNRFromNumRecordsGauss=function(n,numRecords){
  #   f=function(x) numRecords-sqrt(4*n/pi)-x*sqrt(2)/pi*(n*atan(sqrt(n))-sqrt(n))-x*16/(sqrt(3)*pi^1.5)*sqrt(n)*(sqrt(1+1/n)*asinh(sqrt(n))-1)
  #   return(uniroot(f,interval=c(-3,3))$root)
  
  SNR=(numRecords-sqrt(4*n/pi))/(sqrt(2)/pi*(n*atan(sqrt(n))-sqrt(n)))#+16/(sqrt(3)*pi^1.5)*sqrt(n)*(sqrt(1+1/n)*asinh(sqrt(n))-1))
  return(SNR)
}  

meanSkewedDistr=function(delta,dtype="Student",nu=3){
  myscale=(delta^2-1)/(delta)
  if(dtype=="Gauss"){
    myscale=myscale*sqrt(2/pi)*2
  }else if(dtype=="Student"){
    myscale=myscale*2*sqrt(nu)*gamma((nu+1)/2)/(sqrt(pi)*(nu-1)*gamma(nu/2))
  }
  return(myscale)  
}

mean2SkewedDistr=function(delta,dtype="Student",nu=3){
  myscale=(1+delta^6)/(delta^2*(1+delta^2))
  if(dtype=="Gauss"){
    myscale=myscale*sqrt(4/pi)*Gamma(3/2)
  }else if(dtype=="Student"){
    myscale=myscale*nu*gamma((nu-2)/2)/(gamma(nu/2))/2
  }
  return(myscale)  
}

rsktNormalized=function(n,nu,delta){
  return((rskt(n,nu,delta)-meanSkewedDistr(delta,dtype="Student",nu=nu))/(sqrt(mean2SkewedDistr(delta,dtype="Student",nu=nu)-meanSkewedDistr(delta,dtype="Student",nu=nu)^2)))
}

toList <- function(x){
  return(as.list(setNames(x,x)))
}
