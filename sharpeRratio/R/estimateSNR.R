library("ghyp")
library("Rcpp")


a=get(load("data/a_vs_R0dN.rda"))

theta=function(R0dN,nu){
  myEstimate=a(R0dN)*(1-8/3*nu^(-1.5))
  return(myEstimate)
}

estimateSNR=function(x,numPerm=1000){

  if(length(x)<1){
    return(NULL)
  }
  x=x[is.finite(x)]
  if(length(x)<1){
    return(NULL)
  }
  
  myfit=tryCatch(fit.tuv(as.numeric(x),silent=TRUE,nu=6),error=function(e){print(e);return(NA)})
  
  if(!class(myfit)=="mle.ghyp"){
    list(nu=NA,SNR=NA,R0bar=NA,N=length(x))
  }
  
  nu=coef(myfit)$nu

  R0bar=computeR0bar(x,numPerm = numPerm)
  if(abs(R0bar)<1){
    R0bar=0
    return(list(nu=nu,SNR=0,R0bar=0,N=length(x)))
  }

  N=length(x)

  SNR_raw=theta(abs(R0bar)/length(x),nu)    #splines have been calibrated for positive SNRs
  SNR_raw=SNR_raw*  (sign(SNR_raw)>0)                   #if(spline(abs(x)))<0, then set SNR to 0 (this occurs when SNR is very small)     
  SNR=SNR_raw*sign(R0bar)                               #then apply sign of R0bar 

  
  
  return(list(nu=nu,SNR=SNR,R0bar=R0bar,N=length(x)))
}

